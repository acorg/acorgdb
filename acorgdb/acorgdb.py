from collections import abc
from functools import reduce
from reprlib import Repr
from string import ascii_uppercase
from typing import Generator, Iterable, Optional, Callable, Union
import csv
import json
import logging
import os
import random
import re


from airium import Airium
import pandas as pd
import yaml


aRepr = Repr()
aRepr.maxstring = 100
repr = aRepr.repr


class CantGenerateSequenceError(Exception):
    pass


class MissingSequenceError(Exception):
    pass


class EmptySequenceError(Exception):
    """An exception when a sequence is unexpectedly empty."""


class SubstitutionFormatError(Exception):
    pass


class MixedPopulationSubstitutionError(Exception):
    pass


class MissingRecordError(Exception):
    pass


class FrozenJSON:
    """
    Read-only facade for navigating a JSON-like object using attribute notation.

    Notes:
        Based on example 19-5 from Fluent Python, Rahmalho (O'Reilly).
    """

    def __init__(self, mapping):
        self._data = dict(mapping)

    def __repr__(self):
        return "FrozenJSON({})".format(repr(self._data))

    def __getattr__(self, name):
        if hasattr(self._data, name):  # handles calls to .keys etc...
            return getattr(self._data, name)
        else:
            try:
                return FrozenJSON.build(self._data[name])
            except KeyError:
                raise AttributeError(name)

    def __dir__(self):
        return list(self._data.keys())

    @classmethod
    def build(cls, obj):
        if isinstance(obj, abc.Mapping):
            return cls(obj)
        elif isinstance(obj, abc.MutableSequence):
            return [cls.build(item) for item in obj]
        else:
            return obj


class Record(FrozenJSON):
    """
    A Record could be anything in the database, e.g. an Antigen, a Serum, an
    Experiment, a Results table.
    """

    _instances = {}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Register instances for Antigens / Sera
        if isinstance(self, (Antigen, Serum)):
            if self.id in self._instances and self != self._instances[self.id]:
                raise ValueError("Different record already exists")
            else:
                self._instances[self.id] = self

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return hash(self.id)

    def __repr__(self):
        return "{self.__class__.__name__}({self._data})"

    def __str__(self):
        return yaml.dump({self.__class__.__name__: self._data}, sort_keys=False)

    def __getattr__(self, name):
        return super().__getattr__(name)


class Antigen(Record):
    "The Antigen class."

    def __repr__(self):
        return f"Antigen({self._data})"

    @property
    def parent(self):
        # Return None if the antigen doesn't have a parent_id
        try:
            self.parent_id
        except AttributeError:
            return None

        # Return the parent instance
        try:
            return self._instances[self.parent_id]
        except KeyError:
            raise MissingRecordError(f"{self.id} missing from Record._instances")

    @property
    def cluster(self):
        try:
            return self.meta.cluster.name
        except AttributeError:
            return None

    @property
    def _ancestor_generator(self) -> Generator["Antigen", None, None]:
        ancestor = self.parent
        while ancestor is not None:
            yield ancestor
            ancestor = ancestor.parent

    @property
    def _children_generator(self) -> Generator["Antigen", None, None]:
        for record in self._instances.values():
            try:
                if record.parent_id == self.id:
                    yield record
            except AttributeError:
                pass

    @property
    def children(self) -> list["Antigen"]:
        """
        1st generation children of this Antigen.
        """
        return list(self._children_generator)

    @property
    def ancestors(self) -> list["Antigen"]:
        """
        This antigens parent, grandparent etc...
        """
        return list(self._ancestor_generator)

    def has_parent_with_seq(self, gene_name: str = "HA") -> bool:
        """
        Does this antigen have a parent with a sequence for a gene?
        """
        try:
            self.parent.sequence(gene_name)
        except (TypeError, AttributeError, ValueError):
            return False
        else:
            return True

    def alt_parent_id(self, gene_name: str = "HA") -> Union[str, None]:
        """
        The parent_id specified in this antigens alterations field for a gene.
        """
        try:
            return next(
                alt for alt in self.alterations if alt.gene == gene_name.upper()
            ).parent_id
        except (TypeError, AttributeError, StopIteration):
            return None

    def alt_parent(self, gene_name: str = "HA") -> Union["Antigen", None]:
        """
        Lookup the parent antigen specified in the alterations field, if any.
        """
        if alteration_parent_id := self.alt_parent_id(gene_name):
            return self._instances[alteration_parent_id]

    def has_alt_parent_with_seq(self, gene_name: str = "HA") -> bool:
        """
        Does this antigen specify a parent_id in it's alteration field, and does that
        alteration parent have a sequence?
        """
        try:
            self.alt_parent(gene_name).sequence(gene_name)
        except (TypeError, AttributeError):
            return False
        else:
            return True

    def sequence(self, gene: str = "HA") -> str:
        """
        The antigen's sequence. If the antigen does not directly have it's own sequence,
        then one is generated form parent sequences plus any substitutions.

        Args:
            gene_name (str): The name of the gene_name. Defaults to "HA".

        Returns:
            str: The sequence of the gene_name.
        """

        gene = gene.upper()

        own_sequence = self._own_sequence(gene)
        subs = self.substitutions(gene)

        # Variables to make the following logic easier to follow (even if converting to
        # bools isn't strictly necessary)
        has_own_seq = bool(own_sequence)
        has_subs = bool(subs)

        if has_own_seq:

            if has_subs:

                if sequence_consistent_with_aa1(own_sequence, subs):

                    logging.debug(
                        f"Returning {self.id}'s own {gene} sequence which is consistent with its substitutions"
                    )

                    return own_sequence

                else:
                    try:
                        return mutate(own_sequence, subs)

                    except ValueError as err:
                        raise ValueError(
                            f"{self.id} sequence inconsistent with all amino acids "
                            f"gained in {subs} and {err}"
                        )

            else:
                return own_sequence

        else:  # antigen doesn't have its own sequence

            parent = self.parent
            has_parent_with_seq = self.has_parent_with_seq(gene)
            alt_parent = self.alt_parent(gene)
            has_alt_parent_with_seq = self.has_alt_parent_with_seq(gene)

            # If an antigen specifies a parent_id in its alterations field then throw an
            # error if that parent doesn't have a sequence.
            if alt_parent and not has_alt_parent_with_seq:

                raise MissingSequenceError(
                    f"{self.id} specifies {alt_parent.id} as a parent in its "
                    f"alterations field, but {alt_parent.id} doesn't have a sequence "
                    "(or one cannot be constructed)\n"
                    f"{alt_parent}"
                )

            elif has_alt_parent_with_seq and has_subs:

                logging.debug(
                    f"{self.id} {gene} sequence generated by mutating alt parent ({alt_parent.id}) "
                    f"with: {subs}"
                )

                return mutate(alt_parent.sequence(gene), subs)

            elif has_parent_with_seq and has_subs:

                logging.debug(
                    f"{self.id} {gene} sequence generated by mutating parent ({parent.id}) with: "
                    f"{subs}"
                )

                return mutate(parent.sequence(gene), subs)

            elif not (has_alt_parent_with_seq or has_parent_with_seq):

                raise ValueError(f"{self.id} doesn't have a parent with a sequence")

            # If the antigen doesn't have any substitutions, but has a parent
            # (or alteration parent) with a sequence, then return the parent sequence
            elif not has_subs and (has_alt_parent_with_seq or has_parent_with_seq):

                if has_alt_parent_with_seq:

                    logging.debug(
                        f"{self.id} has no substitutions. It's {gene} sequence was taken directly "
                        f"from that of its alt_parent: {alt_parent.id}"
                    )

                    return alt_parent.sequence(gene)

                else:

                    logging.debug(
                        f"{self.id} has no substitutions. It has no {gene} alt parent. It's "
                        f"{gene} sequence was taken directly its parent: {parent.id}"
                    )

                    return parent.sequence(gene)

            else:
                raise CantGenerateSequenceError()

    def substitutions(self, gene_name: str) -> Union[list[str], None]:
        """Substitutions associated directly with this antigen."""
        try:
            for alteration in self.alterations:
                if alteration.gene.upper() == gene_name.upper():
                    return alteration.substitutions
        except AttributeError:
            return None

    def _own_sequence(self, gene_name: str) -> Union[str, None]:
        """A sequence directly associated with this record."""
        try:
            for gene in self.genes:
                if gene.gene.upper() == gene_name.upper():
                    return gene.sequence
        except AttributeError:
            return None

    def _repr_html_(self, seen=None):
        """
        HTML view of the antigen with expandable hierarchy.

        Args:
            seen (set): A set of ids that have already been seen. This is used to
                prevent infinite loops if antigens have circular ancestries.
        """
        seen = set() if seen is None else seen

        seen.add(self.id)

        a = Airium()

        with a.h4():
            a(f"Antigen: {self.id}")

        with a.ul(style="list-style-type:none;"):

            for k, v in self._data.items():

                # Handle alterations separately incase they include a parent_id that
                # should be expandable
                if k == "alterations":
                    with a.li():
                        a("<strong>alterations:</strong>")

                    with a.ul(style="list-style-type:none;"):
                        for alteration in v:

                            with a.li():

                                if "parent_id" in alteration:
                                    alt_id = alteration["parent_id"]

                                    # Keep track of which alt_ids have been seen
                                    # so we avoid infinite loops
                                    if alt_id not in seen:
                                        alt_parent = self._instances[alt_id]
                                        expandable_html(
                                            label=repr(alteration),
                                            content=alt_parent._repr_html_(seen=seen),
                                            a=a,
                                        )
                                        seen.add(alt_id)

                                    else:
                                        with a.p():
                                            a((repr(alteration)))
                                        with a.p():
                                            a(
                                                f"<strong>WARNING loop created</strong> "
                                                f"({alt_id} seen before)"
                                            )

                                else:
                                    a(repr(alteration))

                # id is shown at top, parent_id is shown at bottom
                elif k not in {"id", "parent_id"}:
                    with a.li():
                        a(f"<strong>{k}:</strong> {repr(v)}")

            if "parent_id" in self._data:

                if self.parent_id not in seen:

                    with a.li():
                        a(f"<strong>parent_id:</strong> {self.parent_id}")

                    expandable_html(
                        label="", content=self.parent._repr_html_(seen=seen), a=a
                    )

                    seen.add(self.parent_id)

                else:
                    with a.p():
                        a(
                            f"<strong>WARNING loop created</strong> ({self.parent_id} seen before)"
                        )

        return str(a)


def expandable_html(label: str, content: str, a: Airium) -> str:
    with a.details(style="margin-left:2em;"):
        with a.summary():
            a(label)
        with a.p():
            a(content)


class Serum(Record):
    "The Serum class."

    @property
    def antigen(self) -> Antigen:
        """The antigen used to raise this serum."""
        return self._instances[self.strain_id]

    def _repr_html_(self):
        """
        HTML view of the serum with expandable hierarchy.
        """
        a = Airium()

        with a.h4():
            a(f"Serum: {self.id}")

        with a.ul(style="list-style-type:none;"):

            for k, v in self._data.items():
                if k not in {"id", "strain_id"}:
                    with a.li():
                        a(f"<strong>{k}:</strong> {repr(v)}<br>")

            if "strain_id" in self._data:
                expandable_html(
                    label=f"<strong>strain_id:</strong> {self.antigen.id}",
                    content=self.antigen._repr_html_(),
                    a=a,
                )

        return str(a)


class Result(Record):
    """The results class."""

    def __repr__(self):
        return """Results(
            antigens ({}): {}
            sera ({}): {}
        )""".format(
            len(self.antigen_ids), self.antigen_ids, len(self.serum_ids), self.serum_ids
        )

    def __getattr__(self, name):

        name = "titers" if name == "titers_wide" else name

        value = super().__getattr__(name)

        if name == "titers":
            df = pd.DataFrame(value, self.antigen_ids, self.serum_ids)
            df.columns.name = "serum"
            df.index.name = "antigen"
            return df

        else:
            return value

    @property
    def titers_long(self):
        df = (
            pd.melt(
                self.titers, ignore_index=False, var_name="serum", value_name="titer"
            )
            .eval(f"file = '{self.file}'")
            .reset_index()
        )

        has_comma_mask = df["titer"].str.contains(",")

        if not has_comma_mask.any():
            return df

        else:
            df_with_comma = df[has_comma_mask].copy()
            df_without_comma = df[~has_comma_mask].copy()

            # .split makes comma delimited elements in titer into lists .explode
            # puts each element in the lists into a new row (duplicating values
            # in other columns)
            df_with_comma["titer"] = df_with_comma["titer"].str.split(",")
            df_with_comma = df_with_comma.explode("titer")

            return pd.concat([df_without_comma, df_with_comma])


class Experiment(Record):
    """The experiment class."""

    def __repr__(self):
        return """Experiment(
            id={},
            name={},
            description={}
        )
        """.format(
            self.id, self.name, self.description
        )

    def __getattr__(self, name):
        if name == "results":
            return [Result(result) for result in self._data[name]]
        else:
            return super().__getattr__(name)

    @property
    def titers_wide(self):
        """
        Concatenate all titer tables, joining tables top to bottom.

        Columns are sera, rows are antigens.
        """
        return pd.concat(result.titers_wide for result in self.results)

    @property
    def titers_long(self) -> pd.DataFrame:
        """
        Concatenate titer tables in long format. Tables are joined top to bottom.
        """
        return (
            pd.concat(result.titers_long for result in self.results)
            .eval(f"experiment_id = '{self.id}'")
            .reset_index(drop=True)
        )


class Database:
    def __init__(
        self,
        antigens: Iterable[Antigen],
        sera: Iterable[Serum],
        experiments: Iterable[Experiment],
    ):

        self.antigens = antigens
        self.sera = sera
        self.experiments = experiments
        self._by_id = {}
        for item in *antigens, *sera, *experiments:
            self._by_id[item.id] = item

    def __repr__(self):
        return "DataBase({}, {}, {})".format(
            repr(self.antigens), repr(self.sera), repr(self.experiments)
        )

    def __getitem__(self, item):
        return self._by_id[item]

    @classmethod
    def from_dir(cls, directory: str) -> "Database":
        """
        Load a database from any directory.

        Args:
            directory (str): The directory should contain 'antigens.json',
                'sera.json' and 'results.json' files in the same format as the
                acorgdb repo.

        Returns:
            Database
        """
        with open(os.path.join(directory, "antigens.json")) as f:
            antigens = tuple(Antigen(a) for a in json.load(f))

        with open(os.path.join(directory, "sera.json")) as f:
            sera = tuple(Serum(a) for a in json.load(f))

        with open(os.path.join(directory, "results.json")) as f:
            experiments = tuple(Experiment(a) for a in json.load(f))

        return cls(antigens=antigens, sera=sera, experiments=experiments)

    @property
    def titers_long(self) -> pd.DataFrame:
        """
        Concatenate titer tables in long format.
        """
        return pd.concat(
            experiment.titers_long for experiment in self.experiments
        ).reset_index(drop=True)

    def add_attr_to_df(self, column_name: str, attr: str, id_column) -> Callable:
        """
        Makes a function for populating DataFrames with record attributes.

        Example:

        >>> df
          antigen   serum titer
        0  AVW18R  854SXY    40
        1  AVW18R  YXR994    80

        >>> df.pipe(db.lookup_record_attr("sr_long", attr="long", id_column="serum"))
          antigen   serum titer                             sr_long
        0  AVW18R  854SXY    40            A/GOOSE/GUANGDONG/1/1996
        1  AVW18R  YXR994    80        A/VIETNAM/1194/2004-NIBRG-14

        Args:
            column_name: Name of the resulting column.
            attr: Name of the attribute to lookup.
            id_column: Name of column in passed DataFrame containing ids.
        """

        def fun(df: pd.DataFrame) -> pd.DataFrame:
            f"Adds {column_name} to a DataFrame using {attr} attribute of {id_column} IDs."
            df[column_name] = [getattr(self[id], attr) for id in df[id_column]]
            return df

        return fun


def generate_id(length=6):
    """
    Generate a database ID.

    Args:
        length (int): Length of the id. Default=6.

    Returns:
        str
    """
    chars = tuple(list(range(10)) + list(ascii_uppercase))
    return "".join(str(random.choice(chars)) for _ in range(length))


def print_csv_as_json(path):
    """
    Print data in a CSV file as JSON.

    Args:
        path (str): Path to CSV file.
    """
    with open(path, "r") as f:
        data = [line for line in csv.reader(f)]
    print(json.dumps(data, indent=4))


def load_jsons(json_paths: Iterable[str], cls: Optional[type] = None) -> dict:
    """
    Load multiple JSON files, returning a single dict with "id" as the key.

    Args:
        json_paths (Iterable[str]): A collection of file paths to JSON files.
        cls (optional): The class to use for creating objects from the JSON records. (E.g. Antigen, Serum).

    Returns:
        dict: A dictionary with "id" as the key and the corresponding JSON record as the value.

    Raises:
        ValueError: If multiple records with the same "id" exist.

    Example:
        >>> json_paths = ['/path/to/file1.json', '/path/to/file2.json']
        >>> data = load_jsons(json_paths, cls=Antigen)
    """
    d = {}
    for path in json_paths:
        with open(path) as fobj:
            for record in json.load(fobj):

                id_ = record["id"]

                if id_ in d:
                    raise ValueError(
                        f"Different records with the same id: {id_}\n"
                        f"Value in {path}:\n{record}\n"
                        f"Other value:\n{d[id_]}\n"
                    )

                d[id_] = record if cls is None else cls(record)

    return d


def load_from_dir(directory: str) -> Database:
    raise NotImplementedError("this is now Database.from_dir")


def load_tables(directory: str, index) -> list[pd.DataFrame]:
    """
    Load tables from the database.

    Args:
        path (str): Passed to load_from_dir.
        index (int): Top level index into results.json to access.

    Returns:
        list containing pd.DataFrame
    """
    _, _, results = load_from_dir(directory)
    return [
        pd.DataFrame(table.titers, table.antigen_ids, table.serum_ids)
        for table in results[index].results
    ]


def mutate(sequence: str, substitutions: list[str]) -> str:
    """
    Mutates a given sequence based on a list of substitutions.

    Args:
        sequence (str): The input sequence to be mutated.
        substitutions (list[str]): A list of substitutions in the format "XnY",
            where X is the character to be replaced, n is the site index (1-based),
            and Y is the character to replace with.

    Returns:
        str: The mutated sequence.

    Raises:
        ValueError: If an unrecognised substitution is encountered or if the
            sequence is inconsistent with a substitution.

    """
    if not substitutions:
        return sequence

    if not sequence:
        raise EmptySequenceError()

    sequence = list(sequence)

    for substitution in substitutions:

        aa0, site, aa1 = substitution_components(substitution)
        index = site - 1

        # Sequence has the old amino acid, so update it
        if sequence[index] == aa0:
            sequence[index] = aa1

        else:
            raise ValueError(
                f"Sequence inconsistent with {substitution}. Expected '{aa0}' but find "
                f"'{sequence[index]}' at {site}."
            )

    return "".join(sequence)


def sequence_consistent_with_aa1(sequence: str, substitutions: list[str]) -> bool:
    """
    Check if a sequence is consistent with all the amino acids gained in substitutions.
    If the amino acid gained in any substitution is inconsistent with the sequence, False
    is returned.

    For example, the sequence "MKTLGD" is consistent with all the amino acids gained
    (aa1s) in the substitutions: ["D1M", "L3T"].

    Args:
        sequence (str): The input sequence to check.
        substitutions (list[str]): The list of substitutions in the format 'XnY'.

    Returns:
        bool: True if the sequence is consistent with the amino acids gained in
            substitutions, False otherwise.
    """
    for substitution in substitutions:
        _, site, char_gained = substitution_components(substitution)
        index = site - 1
        if sequence[index] != char_gained:
            return False
    else:
        return True


def substitution_components(substitution: str) -> tuple[str, int, str]:
    """
    Extracts the components from a substitution string.

    Args:
        substitution (str): The substitution string to extract components from.

    Returns:
        tuple[str, int, str]: A tuple containing the extracted components:
                              - The first character of the substitution string.
                              - The numeric portion of the substitution string.
                              - The last character of the substitution string.

    Raises:
        MixedPopulationSubstitutionError: If the substitution contains a mixed population
            amino acid (e.g. "N145N-K").
        SubstitutionFormatError: If the substitution string is not in the expected
            format.
    """
    if match := re.match(r"^(\w)(\d+)(\w)$", substitution):
        aa0, site, aa1 = match.groups()
        if aa0 == aa1:
            raise SubstitutionFormatError(f"aa lost matches aa gained: {substitution}")
        return aa0, int(site), aa1
    else:
        if (
            re.match(r"^\w\d+\w-\w$", substitution)  # aa1 is mixed
            or re.match(r"^\w-\w\d+\w$", substitution)  # aa0 is mixed
            or re.match(r"^\w-\w\d+\w-\w$", substitution)  # both mixed
        ):
            raise MixedPopulationSubstitutionError(substitution)
        else:
            raise SubstitutionFormatError(substitution)


def get_own_subs(ag: Antigen) -> set[str]:
    """
    Get substitutions directly associated with this antigen.

    Args:
        ag: acorgdb.Antigen
    """
    subs = ag.substitutions("HA")
    return set(subs) if subs is not None else set()


def get_ancestor_subs(ag: Antigen) -> set[str]:
    """
    Get substitutions from this antigen's ancestors.

    Args:
        ag: acorgdb.Antigen
    """
    subs = [get_own_subs(anc) for anc in ag.ancestors]
    return set(reduce(set.union, subs)) if subs else set()


def get_sub_pos(sub: str) -> int:
    """
    Get a substitution position.
    """
    match = re.match(r"^[A-Z](\d+)[A-Z](?:\-[A-Z])?$", sub)
    return int(match.groups()[0])


def get_subs_in_name(ag: Antigen | str) -> set[str]:
    """
    Find substitutions in an antigen's long name:

    Args:
        ag: An acorgdb.Antigen or string.
    """
    pattern = r"(?<=-|\/|_| )([A-Z]\d+[A-Z](?:-[A-Z])?)(?=-|\/|_| |$)"
    string = ag.long if isinstance(ag, Antigen) else ag
    return set(re.findall(pattern, string))


def remove_mixed_subs(subs: set[str]) -> set[str]:
    """
    Substitutions from mixed virus populations look like "XnY-Z". This function removes
    any substitution that matches this format.
    """
    return set(sub for sub in subs if not re.match(r"[A-Z]\d+[A-Z]\-[A-Z]$", sub))
