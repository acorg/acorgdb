import regex as re
import unittest

import acorgdb as adb
from acorgdb import Antigen


class TestMutate(unittest.TestCase):
    def test_mutate_single_substitution(self):
        self.assertEqual("NTTRG", adb.mutate("NKTRG", ["K2T"]))

    def test_mutate_multiple_substitutions(self):
        self.assertEqual("NTTRP", adb.mutate("NKTRG", ["K2T", "G5P"]))

    def test_mutate_no_substitutions(self):
        self.assertEqual("NKTRG", adb.mutate("NKTRG", []))

    def test_mutate_substitutions_none(self):
        self.assertEqual("NKTRG", adb.mutate("NKTRG", None))

    def test_mutate_empty_sequence(self):
        with self.assertRaises(adb.EmptySequenceError):
            adb.mutate("", ["K2T", "G5P"])

    def test_inconsistent_substitutions(self):
        with self.assertRaisesRegex(ValueError, "sequence inconsistent with K3P"):
            adb.mutate("ACTGN", ["K3P"])


class TestAntigenSequence(unittest.TestCase):

    def setUp(self) -> None:
        """
        New _instances dictionary so that tests are independent of each other.
        """
        Antigen._instances = dict()

    def test_own_sequence(self):
        """
        Antigen has a sequence and no substitutions (case 1).
        """
        ag = Antigen(
            {
                "id": "ag",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVDTIME"}],
                "wildtype": False,
            }
        )
        self.assertEqual("DQICIGYHANNSTEQVDTIME", ag.sequence("HA"))

    def test_own_sequence_plus_subs(self):
        """
        Antigen has a sequence and substitutions (case 2).
        """
        ag = Antigen(
            {
                "id": "ag",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVDTIME"}],
                "alterations": [
                    {"gene": "HA", "substitutions": ["D1K", "G6T"]},
                    {"gene": "NA", "substitutions": ["D1T", "G6D"]},
                ],
                "wildtype": False,
            }
        )
        self.assertEqual("KQICITYHANNSTEQVDTIME", ag.sequence("HA"))

    def test_parent_has_seq(self):
        """
        Antigen has substitutions, and the parent has a sequence (case 3).
        """
        Antigen(
            {
                "id": "parent",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVQTIME"}],
                "wildtype": False,
            }
        )
        ag = Antigen(
            {
                "id": "child",
                "alterations": [
                    {"gene": "HA", "substitutions": ["D1K", "G6T"]},
                ],
                "parent_id": "parent",
                "wildtype": False,
            }
        )
        self.assertEqual("KQICITYHANNSTEQVQTIME", ag.sequence("HA"))

    def test_antigen_has_parent_no_seq_or_subs(self):
        """
        Antigen with a parent but without a sequence or subs is not implemented (case 4).
        """
        Antigen(
            {
                "id": "parent",
                "wildtype": False,
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVQTIME"}],
            }
        )
        ag = Antigen({"id": "CHILD5", "parent_id": "parent", "wildtype": False})
        msg = (
            "Generating a sequence for an antigen with a parent but without "
            "substitutions should return the parent sequence."
        )
        self.assertTrue(ag.has_parent_with_seq("HA"))
        self.assertFalse(ag.has_alt_parent_with_seq("HA"))
        self.assertEqual(ag.sequence("HA"), ag.parent.sequence("HA"), msg)

    def test_parent_and_antigen_have_sequence(self):
        """
        Where both the parent and antigen have a sequence, the antigen's sequence should
        be used (case 5).
        """
        Antigen(
            {
                "id": "parent",
                "genes": [{"gene": "HA", "sequence": "STEQVQTIME"}],
                "wildtype": False,
            }
        )
        ag = Antigen(
            {
                "id": "child",
                "parent_id": "parent",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHAN"}],
                "wildtype": False,
            }
        )
        self.assertEqual("DQICIGYHAN", ag.sequence("HA"))

    def test_antigen_has_seq_and_subs_and_parent(self):
        """
        Antigen has a sequence, substitutions and a parent (case 6).
        """
        Antigen(
            {
                "id": "parent",
                "genes": [{"gene": "HA", "sequence": "XXX"}],
                "wildtype": False,
            }
        )
        ag = Antigen(
            {
                "id": "child",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVQTIME"}],
                "alterations": [{"gene": "HA", "substitutions": ["D1K", "G6T"]}],
                "parent_id": "parent",
                "wildtype": False,
            }
        )
        self.assertEqual("KQICITYHANNSTEQVQTIME", ag.sequence("HA"))

    def test_antigen_without_parent_or_sequence(self):
        """
        An antigen without a parent or a sequence should raise an error (case 7).
        """
        ag = Antigen({"id": "child", "wildtype": False})
        msg = "child doesn't have a parent"
        with self.assertRaisesRegex(ValueError, msg):
            ag.sequence("HA")

    def test_grandparent_has_sequence(self):
        """
        Case where a grandparent has a sequence, and children have alterations.

        Amino acid at 3 gets altered first to M and then to L.
        """
        Antigen(
            {
                "id": "grandparent",
                "genes": [
                    {"gene": "HA", "sequence": "WSYIVEKINPANDLCYPGNFNDYEELKHLLSR"}
                ],
                "wildtype": False,
            }
        )
        Antigen(
            {
                "id": "parent",
                "parent_id": "grandparent",
                "alterations": [{"gene": "HA", "substitutions": ["Y3M", "N19T"]}],
                "wildtype": False,
            }
        )
        ag = Antigen(
            {
                "id": "CHILD3",
                "parent_id": "parent",
                "alterations": [{"gene": "HA", "substitutions": ["M3L"]}],
                "wildtype": False,
            }
        )
        self.assertEqual("WSLIVEKINPANDLCYPGTFNDYEELKHLLSR", ag.sequence("HA"))

    def test_antigen_that_specifies_aa1s_present(self):
        """
        Antigen lists substitutions and a sequence. All the substitutions and the amino
        acids that are gained in these substitutions are already present in it's
        sequence.
        """
        ag = Antigen(
            {
                "id": "child",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVQTIME"}],
                "alterations": [
                    {"gene": "HA", "substitutions": ["K1D", "T6G", "D21E"]}
                ],
                "wildtype": False,
            }
        )
        self.assertEqual("DQICIGYHANNSTEQVQTIME", ag.sequence("HA"))

    def test_antigen_specifies_inconsistent_substitution(self):
        """
        Like above, but the sequence has an E at 21 and the substitution at site 21
        gains a K. (Amino acids gained in other substitutions all match the sequence).
        If not all substitution aa1s are consistent with the sequence, a ValueError
        should be raised.
        """
        ag = Antigen(
            {
                "id": "child",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVQTIME"}],
                "alterations": [
                    {"gene": "HA", "substitutions": ["K1D", "T6G", "D21K"]}
                ],
                "wildtype": False,
            }
        )
        msg = (
            "child sequence inconsistent with all amino acids gained in "
            r"\['K1D', 'T6G', 'D21K'\] and sequence inconsistent with K1D"
        )

        with self.assertRaisesRegex(ValueError, msg):
            ag.sequence("HA")

    def test_passing_unrecognised_substitution(self):
        """
        Passing an unrecognised substitution should raise an appropriate error.
        """
        Antigen(
            {
                "id": "parent",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVDTIME"}],
                "wildtype": False,
            }
        )
        ag = Antigen(
            {
                "id": "child",
                "parent_id": "parent",
                "alterations": [
                    {"gene": "HA", "substitutions": ["D1D-K"]},
                ],
                "wildtype": False,
            }
        )
        with self.assertRaises(adb.MixedPopulationSubstitutionError):
            ag.sequence("HA")


class TestAntigenSequenceParentSpecifiedInAlterations(unittest.TestCase):
    """
    In some cases an antigen's parent may not have a sequence and the alterations field
    might contain the ID of the parent sequence.

    For example GKHDPB's parent (WU29LG) does not have a sequence. The HA and NA
    sequences are stored in the alterations field:

    "alterations": [
              {
                  "gene": "HA",
                  "parent_id": "II6SL4"
              },
              {
                  "gene": "NA",
                  "parent_id": "II6SL4"
              }
          ],

    These tests check this functionality.
    """

    def setUp(self) -> None:
        """
        New _instances dictionary so that tests are independent of each other.
        """
        Antigen._instances = dict()

    def test_uses_alteration_parent_by_default(self):
        """
        If a parent_id is present in the alterations field, this parent should be used as
        the sequence source by default.
        """
        Antigen(
            {
                "id": "parent1",
                "genes": [{"gene": "HA", "sequence": "SHOULDBEAVOIDED"}],
            }
        )
        Antigen(
            {
                "id": "parent2",
                "genes": [{"gene": "HA", "sequence": "KEQUENCETOUSE"}],
            }
        )
        # Construct a child that specifies a main parent AND a parent in it's alterations
        child = Antigen(
            {
                "id": "child1",
                "parent_id": "parent1",
                "alterations": [
                    {
                        "gene": "HA",
                        "substitutions": ["K1S"],
                        "parent_id": "parent2",
                    },
                ],
            }
        )
        self.assertEqual(
            "SEQUENCETOUSE",
            child.sequence("HA"),
            "child antigen's sequence is being constructed from the parent_id specified "
            "in the alterations field, but should be taken from the antigen's main "
            "parent",
        )

    def test_falls_back_to_main_parent(self):
        """
        If the alterations field lacks a parent_id then use the main parent id.
        """
        Antigen({"id": "parent", "genes": [{"gene": "HA", "sequence": "QRSTUVWXYZ"}]})

        child = Antigen(
            {
                "id": "child",
                "parent_id": "parent",
                "alterations": [{"gene": "HA", "substitutions": ["Q1K"]}],
            }
        )

        self.assertEqual("KRSTUVWXYZ", child.sequence("HA"))

    def test_alt_parent_specified_but_alt_parent_lacks_sequence(self):
        """
        If an antigen specifies a parent_id in it's alterations field, and that parent
        doesn't have a sequence, an error should be raised (even if a main parent exists
        with a sequence).
        """
        Antigen(
            {"id": "main_parent", "genes": [{"gene": "HA", "sequence": "QRSTUVWXYZ"}]}
        )
        Antigen({"id": "alt_parent"})
        child = Antigen(
            {
                "id": "child",
                "parent_id": "main_parent",
                "alterations": [
                    {
                        "gene": "HA",
                        "substitutions": ["Q1K"],
                        "parent_id": "alt_parent",
                    },
                ],
            }
        )
        with self.assertRaisesRegex(
            ValueError, "alt_parent doesn't have a parent with a sequence"
        ):
            child.sequence("HA")

    def test_alt_parent_id(self):
        """
        Check that a parent_id in the alterations field is looked up correctly.
        """
        ag = Antigen(
            {
                "id": "child",
                "alterations": [
                    {"gene": "HA", "parent_id": "altparent"},
                ],
            }
        )
        self.assertEqual("altparent", ag.alt_parent_id("HA"))


class TestSubstitutionComponents(unittest.TestCase):
    def test_k1d(self):
        self.assertEqual(("K", 1, "D"), adb.substitution_components("K1D"))

    def test_t6g(self):
        self.assertEqual(("T", 6, "G"), adb.substitution_components("T6G"))

    def test_d21e(self):
        self.assertEqual(("D", 21, "E"), adb.substitution_components("D21E"))

    def test_trailing_characters(self):
        with self.assertRaises(adb.MixedPopulationSubstitutionError):
            adb.substitution_components("A45T-I")

    def test_leading_characters(self):
        with self.assertRaises(adb.MixedPopulationSubstitutionError):
            adb.substitution_components("A-A45T")

    def test_amino_acids_differ(self):
        with self.assertRaises(adb.SubstitutionFormatError):
            adb.substitution_components("A12A")


class TestAntigen(unittest.TestCase):
    def test_wildtype_antigen_with_parent(self):
        """
        A wildtype is allowed to have a parent.
        """
        Antigen({"id": "AGTEST1_PARENT", "wildtype": True})
        Antigen(
            {"id": "AGTEST1", "wildtype": True, "parent_id": "AGTEST1_PARENT"}
        ).parent

    def test_missing_parent_instance(self):
        """
        If an antigen has a parent_id that doesn't exist in Records._instances, a
        MissingParentError should be raised.
        """
        ag = Antigen(
            {"id": "AGTEST3", "parent_id": "AGTEST3_PARENT", "wildtype": False}
        )
        with self.assertRaises(adb.MissingRecordError):
            ag.parent

    def test_no_parent_id(self):
        """
        An antigen without a parent_id should have None as it's parent attribute.
        """
        ag = Antigen({"id": "AGTEST3", "wildtype": True})
        self.assertIsNone(ag.parent)


class TestGatesH5Substitutions(unittest.TestCase):
    """
    Tests for the H5 database.

    Would put these somewhere in lib/test/test_db.py. But (a) Sam removed the mechanism
    I'd written for writing bespoke tests for a specific experiment and (b) the following
    tests use the functionality in the acorgdb module I've written.

    At the moment I'm just testing the Gates antigens.
    """

    @classmethod
    def setUpClass(cls) -> None:
        db = adb.load(("h5", "experiments", "h5_mutants"))
        cls.antigens = [ag for ag in db.antigens if ag.id in GATES_AG_IDS]

    def test_antigen_sequence(self):
        """
        Should be able to generate a sequence for all antigens.
        """
        for ag in self.antigens:
            with self.subTest(ag_id=ag.id):

                try:
                    ag.sequence("HA")
                except adb.MixedPopulationSubstitutionError:
                    # Ignore errors thrown by substitutions involving mixed populations
                    # of amino acids
                    ...

    def test_subs_in_own_alterations_present_in_long_name(self):
        """
        Substitutions in an antigens own alterations field should be in the antigen's
        long name. (Long name may also include additional substitutions from parents.)
        """
        for ag in self.antigens:
            with self.subTest(ag_id=ag.id):

                if ag.id == "IWY9GS":
                    # This is the antigen with long name "NODE2"
                    continue

                own_subs = adb.get_own_subs(ag)
                subs_in_name = adb.get_subs_in_name(ag)

                subs_in_alterations_but_not_name = own_subs - subs_in_name

                if all(
                    re.match(r"^\w\d+\w-\w$", s)
                    for s in subs_in_alterations_but_not_name
                ):
                    # If all substitutions that differ involve mixed amino acids, then
                    # skip the test
                    continue
                else:
                    self.assertFalse(
                        subs_in_alterations_but_not_name,
                        f"\n{ag.id} has substitutions in alterations not in its name\n"
                        f"alterations: {sorted(own_subs, key=adb.get_sub_pos)}\n"
                        f"-      name: {sorted(subs_in_name, key=adb.get_sub_pos)}\n"
                        "----------------------------\n"
                        f"difference : {sorted(subs_in_alterations_but_not_name, key=adb.get_sub_pos)}\n\n"
                        f"{ag.long}",
                    )

    def test_subs_in_long_name_present_in_alterations(self):
        """
        All substitutions in an antigen's name should be in its, and its parents,
        alterations.
        """
        for ag in self.antigens:
            with self.subTest(ag_id=ag.id):

                own_subs = adb.get_own_subs(ag)
                anc_subs = adb.get_ancestor_subs(ag)

                subs_in_name = adb.get_subs_in_name(ag)

                subs_in_name_not_alterations = subs_in_name - (own_subs | anc_subs)

                if all(
                    re.match(r"^\w\d+\w-\w$", s) for s in subs_in_name_not_alterations
                ):
                    # If all substitutions that differ involve mixed amino acids, then
                    # skip the test
                    continue
                else:
                    self.assertFalse(
                        subs_in_name_not_alterations,
                        f"\n{ag.id} has substitutions in its name but not in alterations\n"
                        f"         name: {sorted(subs_in_name,key=adb.get_sub_pos)}\n"
                        f"- alterations: {sorted(own_subs.union(anc_subs), key=adb.get_sub_pos)}\n"
                        "----------------------------\n"
                        f"difference : {sorted(subs_in_name_not_alterations, key=adb.get_sub_pos)}\n\n"
                        f"{ag.long}",
                    )


class TestRemoveMixedSubs(unittest.TestCase):
    def test_remove_mixed_subs(self):
        self.assertEqual(
            {"N145K"},
            adb.acorgdb.remove_mixed_subs({"R189R-G", "K267K-I", "I213I-V", "N145K"}),
        )


class TestGetSubsInName(unittest.TestCase):

    def test_case_a(self):
        name = "NODE2-PR8_A/WHOOPERSWAN/MONGOLIA/244/2005NA-HA-K140R/S155P/R189V"
        expect = {"K140R", "S155P", "R189V"}
        self.assertEqual(expect, adb.get_subs_in_name(name))

    def test_mixed_subs(self):
        name = (
            "NODE2-PR8_A/WHOOPERSWAN/MONGOLIA/244/2005NA-HA-K140K-S/S155S-L/"
            "R189R-W-HA-N87N-Y/I151I-T/A156A-T/N165N-K/N220N-Y/A238A-T"
        )
        expect = {
            "N87N-Y",
            "K140K-S",
            "S155S-L",
            "R189R-W",
            "I151I-T",
            "A156A-T",
            "N165N-K",
            "N220N-Y",
            "A238A-T",
        }
        self.assertEqual(expect, adb.get_subs_in_name(name))


GATES_AG_IDS = {
    "CR426H",
    "4IYM4V",
    "Z0URGP",
    "A60614",
    "RAEMYS",
    "MW8V63",
    "Q849HU",
    "ZCJW2Y",
    "23H9PF",
    "9FJOVM",
    "JXWXB9",
    "42PCG1",
    "B40617",
    "ZEX3QF",
    "O75M60",
    "X0G8RF",
    "0WBVE7",
    "BB5T4D",
    "PN7RPV",
    "E3N9D2",
    "4ACF1F",
    "UZ453S",
    "A020RM",
    "IWY9GS",
    "SMBLR3",
    "G39QDE",
    "RUHSHZ",
    "MNXI7R",
    "SG915J",
    "AL7DFO",
    "2U7GA8",
    "M1YOA4",
    "AVW18R",
    "4OG5AW",
    "KFZMWF",
    "077RCO",
    "S609T5",
    "J6769R",
    "AD80EE",
    "NF4M95",
    "CZII7I",
    "FYGNDL",
    "CDC744",
    "RNIOE9",
    "Q1FOTE",
    "HR7ADG",
    "M3EAKZ",
    "Q5NWOP",
    "O5HVTQ",
    "P64HM1",
    "GSE5V1",
    "XCP0L9",
    "0E1729",
    "WLTTW6",
    "62B056",
    "976UL8",
    "VDBRMJ",
    "ZKS5VL",
    "YUB7LP",
    "PHBEOH",
    "NK4KU7",
    "U5D8XV",
    "ZMOGC6",
    "YAZ9G9",
    "LEOOFC",
    "50JGPD",
    "1ATQHR",
    "GOFCZ2",
    "NVY8TY",
    "XIZOZL",
    "FREONF",
    "2E5FYV",
    "II9N0F",
    "LY1ZXH",
    "LOK77D",
    "YWYTH0",
    "EAR84J",
    "LN270S",
    "YGFBY1",
    "0E8263",
    "JM6CME",
    "863GD2",
    "YRRSQW",
    "SGE32R",
    "LGIN8G",
    "SYPMQQ",
    "9D44Z5",
    "CYADGD",
    "C6NT71",
    "NG5QE4",
    "7RF735",
    "OJUULM",
    "T0019F",
    "DZCN90",
    "JBSRTB",
    "7CJ09K",
    "MSBJ23",
    "7O3ND7",
    "E2C546",
    "BHFIPQ",
    "79PVXL",
    "7R1NOM",
    "J67DSZ",
    "ZGDMT6",
    "IJG4A0",
    "X9MYQG",
    "6BJZ75",
    "98BVFS",
    "Z9FCCV",
    "BVQCLH",
    "CXSYDG",
    "YK7L72",
    "E1PSP7",
    "BYM44A",
    "HLAZAS",
    "BS0GTL",
    "87RFX0",
    "FY85GB",
    "BB7QI4",
    "HRBH9E",
    "5LFZYH",
    "F8XVZS",
    "96MLNJ",
    "SGYS4T",
    "YES8QV",
    "XILW2X",
    "OK7B3I",
    "ZLLK0Z",
    "JVOM0A",
    "QXJ2RR",
    "B66CDA",
    "2F8BUR",
    "0NPYMF",
    "336895",
    "14846I",
    "R9U5U9",
    "0SCMX5",
    "MDIAUC",
    "09JE2G",
    "EVZH8R",
    "N17BI0",
    "ABE133",
    "QSQ2MP",
    "3ZJI4D",
    "F4340C",
    "4LGQXJ",
    "NNKMFK",
    "ALY0C1",
    "I4Z8K3",
    "XRGD0G",
    "JJPPK0",
    "EI12H4",
    "M8CWF8",
    "7CMDDI",
    "J2081A",
    "0SNIT2",
    "9NM7VY",
    "XZERFV",
    "J00PN7",
    "1DZ6GC",
    "VI9XMY",
    "LQPS0R",
    "69E969",
    "EMKZON",
    "FLU02U",
    "DTLHMA",
    "57E5F1",
    "WIH89H",
    "FAZNQK",
    "8XZWB7",
    "EKGFE6",
    "7WT339",
    "IGA2EK",
    "F6BIV9",
    "8VDHV3",
    "ZZWLU2",
    "ROXNKV",
    "PA9ZQW",
    "F4EFF6",
    "HI2QT7",
    "TSL8H3",
    "QTX4G8",
    "7P6Z0F",
    "NX80D9",
    "1A694W",
    "J1B2B3",
    "MT1ORR",
    "EBZJTC",
    "2A02UY",
    "ATNRBU",
    "PXLY8Y",
    "QEQ7VN",
    "YZDE94",
    "Z8CUKN",
    "NVIDW4",
    "SL0DT0",
    "M5HQ8H",
    "N011M3",
    "70UY68",
    "ZEPMD6",
    "CGOCNM",
    "ITAXYO",
    "YNG3R8",
    "FD6L5U",
    "8Z42TK",
    "6PW9HP",
    "QLU477",
    "JVNZ03",
    "3SLPWJ",
    "5LJB91",
    "7ZJT88",
    "ZPQ8B1",
    "G20RNU",
    "F7L9KF",
    "PMQ2YR",
    "ZSUSHW",
    "G36UCH",
    "EO3354",
    "ARTLF7",
    "QESLLA",
    "J1U671",
    "LXITI6",
    "LM9DZU",
    "HNB77P",
    "97E9FA",
    "6PKWVC",
    "PPELIG",
    "V3ZGX2",
    "UXD4C6",
    "KFBWCS",
    "Q3SNN4",
    "AH1K7V",
    "LU0OG0",
    "L5IHLA",
    "T5T5TP",
    "3ZUMYS",
    "LZTZVX",
    "3MF08T",
    "AV15CU",
    "JU4FXO",
    "RENR39",
    "U18MVX",
    "GZRBXD",
    "I7HAMW",
    "7N90YJ",
    "IFWA7P",
    "19X1LI",
    "DA3VRS",
    "46993C",
    "2H9FH2",
    "B32HA5",
    "6A678E",
    "GYR6NP",
    "499788",
    "GJONQU",
    "IRJ359",
    "SHW7E4",
    "9S7QH9",
    "W1MU5U",
    "MFY3MN",
    "6ZXG4P",
    "1E08A5",
    "DIZTCH",
    "DHC1P8",
    "TLI2C6",
    "3TV885",
    "CZQQN8",
    "M2HF0A",
    "UCBNRY",
    "1GP88A",
    "AFJAPC",
    "6A4249",
    "P312BA",
    "C7LJMC",
    "MBBIFN",
    "R9E3TZ",
    "LUU31P",
    "9R2MPL",
    "F5C82S",
    "EAE904",
    "OJUKLX",
    "QR3DFJ",
    "HB3ZJT",
    "Q2TKDB",
    "79JWSM",
    "S31W0J",
    "H7FRDC",
    "J6O6Z1",
    "H9FD8F",
    "Z91HCG",
    "8AFB19",
    "CTDUE4",
    "60Q5IV",
    "LWNOE7",
    "YXRD5E",
    "9MW3CS",
    "QQ5Y1R",
    "PTGY6M",
    "PMBZ5E",
    "PXC55M",
    "YR5LTM",
    "AP243S",
    "3AF041",
    "N6U3UG",
    "G4GF1T",
    "CB2911",
    "F63MRJ",
    "8L0VMP",
    "PTL40K",
    "339238",
    "922D4D",
    "D0B42D",
    "0BA9F6",
    "F7NCLT",
    "06IEHH",
    "BUH31Q",
    "BJ5D7B",
    "LJ5U0L",
    "JA0GNC",
    "P55FV2",
    "P98S0K",
    "110E37",
    "IV3C8Z",
    "81C428",
    "9S3Z11",
    "UAXCC4",
    "04JXFO",
    "I79CNS",
    "G2GJS7",
    "ZBR96W",
    "RZUV0Q",
    "758A38",
    "ULWCIQ",
    "AJJY25",
    "SYU87C",
    "30IP0N",
    "N9JWDE",
    "EH1BNS",
    "7V9VCB",
    "3G73XP",
    "9RSMLI",
    "2P01MY",
    "3029UG",
    "6F73OE",
    "NWYKXW",
    "NA90EM",
    "VLF0Y0",
    "EVSZ1Q",
    "PLQ4E7",
    "L866AW",
    "STHSC9",
    "6LK8SM",
    "ISTT3M",
    "JJGZXB",
    "JEXKOU",
    "OHXJ54",
    "58FF7D",
    "3YWYLU",
    "0J0OQX",
    "87WUQ3",
    "EJE658",
    "DP10XE",
    "582O02",
    "RV2YKX",
    "ZP35P6",
    "NZ7RR9",
    "LBNZS1",
    "7ZDX3D",
    "X1TS4V",
    "9B5HYH",
    "1T9VIZ",
    "F5Q3C4",
    "VC99ZD",
}


if __name__ == "__main__":
    unittest.main()
