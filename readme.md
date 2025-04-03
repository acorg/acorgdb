# acorgdb

Acorg antigen, antisera and titration databases.

This module provides easy access to the database in python.

## Install

Clone and cd into `acorg/acorgdb` then `pip install .` (`pip install -e .[dev]`
if you want to make changes.)

## Example usage

```python
import acorgdb

db = acorgdb.load(("h3n2", "experiments", "SI87_BE92_singles"))

# Access the antigens, sera and experiments associated with a database directory
db.antigens[25]

# You can also access specific antigens, sera and experiments via their ID:
antigen = db["09V8SO"]

# Access attributes of antigens, sera, and experiments:
antigen.isolation.cell

experiment = db["XDGN6Z"]
experiment.name

# Experiments have titers in long and wide format as pandas DataFrames:
db.experiments[0].titers_long

db["XDGN6Z"].titers_wide
```


## Sequences

```python
import acorgdb

db = acorgdb.load(("h5", "experiments", "h5_mutants"))
```

This antigen doesn't have its own sequence:

```python
ag = db["DHC1P8"]
print(ag)
```
```
Antigen:
  id: DHC1P8
  parent_id: IWY9GS
  long: NODE2-PR8_A/WHOOPERSWAN/MONGOLIA/244/2005NA-HA-K140L/S155G/R189I
  wildtype: false
  alterations:
  - gene: HA
    substitutions:
    - K140L
    - S155G
    - R189I
```
but its parent does:
```python
print(ag.parent)
```
```
Antigen:
  id: IWY9GS
  parent_id: TRRDQG
  long: NODE2
  wildtype: false
  alterations:
  - gene: HA
    substitutions:
    - L71I
    - I83A
    - R140K
    ...
  genes:
  - gene: HA
    sequence: DQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDLDGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKANPANDL...

```
The child antigen's sequence is constructed from the parent's sequence while
incorporating the child's substitutions:

```python
print(ag.sequence("HA"))
```
```
DQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDLDGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKANPANDL...
```

If an antigen's parent doesn't have it's own sequence then the parent's parent is
checked etc... until an ancestor is found with a sequence. Substitutions then are
incorporated at each generation until the sequence of interest is generated.

### Sequences with substitutions already incorporated

Sometimes antigens list substitutions that are inconsistent with the parent sequences. 
For example, the substitution might be D1K but the sequence might start PMT... Here, site
1 does not have a D, so there is an inconsistency.

Often, mutants list their substitutions _and_ have a sequence listed. In these cases if
the amino acid _gained_ in a substitution is consistent with sequence position, and this 
is true for all substitutions, then no error is raised. When this is checked, if even a 
single substitution has an amino acid that is gained that is inconsistent with the
sequence, an error is raised. These tests capture this behaviour:

```python
class TestAntigenSequence(unittest.TestCase):
    
    ...
    
    def test_antigen_that_specifies_aa1s_present(self):
        """
        Antigen lists substitutions and a sequence. All the substitutions and the amino
        acids that are gained in these substitutions are already present in it's
        sequence.
        """
        ag = adb.Antigen(
            {
                "id": "CHILD8",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVQTIME"}],
                "alterations": [
                    {"gene": "HA", "substitutions": ["K1D", "T6G", "D21E"]}
                ],
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
        ag = adb.Antigen(
            {
                "id": "CHILD8",
                "genes": [{"gene": "HA", "sequence": "DQICIGYHANNSTEQVQTIME"}],
                "alterations": [
                    {"gene": "HA", "substitutions": ["K1D", "T6G", "D21K"]}
                ],
            }
        )
        msg = (
            "CHILD8 sequence inconsistent with all amino acids gained in "
            r"\['K1D', 'T6G', 'D21K'\] and sequence inconsistent with K1D"
        )

        with self.assertRaisesRegex(ValueError, msg):
            ag.sequence("HA")
```