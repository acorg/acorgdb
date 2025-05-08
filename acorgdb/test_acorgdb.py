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
        with self.assertRaisesRegex(ValueError, "Sequence inconsistent with K3P."):
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
        msg = "child sequence inconsistent with all amino acids gained"

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

    def test_has_alt_parent_with_seq_positive(self):
        """
        Positive test case for Antigen.has_alt_parent_with_seq
        """
        Antigen({"id": "parent"})
        Antigen({"id": "alt_parent", "genes": [{"gene": "HA", "sequence": "QNPS"}]})
        ag = Antigen(
            {
                "id": "child",
                "parent_id": "parent",
                "alterations": [
                    {"gene": "HA", "parent_id": "alt_parent", "substitutions": ["N2K"]}
                ],
            }
        )
        self.assertTrue(ag.has_alt_parent_with_seq("HA"))

    def test_has_alt_parent_with_seq_negative(self):
        """
        Negative test case for Antigen.has_alt_parent_with_seq
        """
        Antigen({"id": "parent"})
        ag = Antigen({"id": "child", "parent_id": "parent"})
        self.assertFalse(ag.has_alt_parent_with_seq("HA"))

    def test_has_parent_with_seq_positive(self):
        """
        Positive test case for Antigen.has_parent_with_seq
        """
        Antigen({"id": "parent", "genes": [{"gene": "HA", "sequence": "QNPS"}]})
        ag = Antigen({"id": "child", "parent_id": "parent"})
        self.assertTrue(ag.has_parent_with_seq("HA"))

    def test_has_parent_with_seq_negative(self):
        """
        Negative test case for Antigen.has_parent_with_seq
        """
        Antigen({"id": "parent"})
        Antigen({"id": "alt_parent", "genes": [{"gene": "HA", "sequence": "QNPS"}]})
        ag = Antigen(
            {
                "id": "child",
                "parent_id": "parent",
                "alterations": [
                    {"gene": "HA", "parent_id": "alt_parent", "substitutions": ["N2K"]}
                ],
            }
        )
        self.assertFalse(ag.has_parent_with_seq("HA"))

    def test_has_parent_without_sequence_but_alt_parent_with_sequence(self):
        """
        In cases where a parent is present (but doesn't have a sequence) and an
        alteration parent is also present (that does have a a sequence), the
        alteration parent should be used for generating a sequence.
        """
        Antigen({"id": "parent"})
        Antigen({"id": "alt_parent", "genes": [{"gene": "HA", "sequence": "QNPS"}]})
        ag = Antigen(
            {
                "id": "child",
                "parent_id": "parent",
                "alterations": [
                    {"gene": "HA", "parent_id": "alt_parent", "substitutions": ["N2K"]}
                ],
            }
        )
        self.assertTrue(ag.has_alt_parent_with_seq("HA"))
        self.assertEqual("QKPS", ag.sequence("HA"))


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


if __name__ == "__main__":
    unittest.main()
