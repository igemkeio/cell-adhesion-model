"""文献値と仮定のメタデータのテスト。

出典のない数値が紛れ込むこと、仮定の記述が欠けることを防ぐ。
"""

import unittest

from adhesion import literature as lit
from adhesion.assumptions import ALL_ASSUMPTIONS, BY_ID, Status, uses


class ReferenceTests(unittest.TestCase):
    def test_every_reference_has_a_doi(self):
        for key, reference in lit.REFERENCES.items():
            with self.subTest(reference=key):
                self.assertTrue(reference.doi)
                self.assertTrue(reference.url.startswith("https://doi.org/"))

    def test_every_reference_has_bibliographic_detail(self):
        for key, reference in lit.REFERENCES.items():
            with self.subTest(reference=key):
                self.assertTrue(reference.authors)
                self.assertTrue(reference.title)
                self.assertTrue(reference.journal)
                self.assertGreater(reference.year, 1900)

    def test_short_citation_format(self):
        self.assertEqual(lit.REFERENCES["costan2024"].short, "Costan et al. (2024)")


class LiteratureValueTests(unittest.TestCase):
    def test_every_value_points_at_a_registered_reference(self):
        for value in lit.ALL_VALUES:
            with self.subTest(symbol=value.symbol):
                self.assertIn(value.source, lit.REFERENCES)

    def test_every_value_describes_its_quantity(self):
        for value in lit.ALL_VALUES:
            with self.subTest(symbol=value.symbol):
                self.assertTrue(value.quantity)
                self.assertTrue(value.symbol)

    def test_unknown_source_is_rejected(self):
        with self.assertRaises(KeyError):
            lit.LiteratureValue(
                symbol="x", value=1.0, unit="m", quantity="test", source="nonexistent"
            )

    def test_values_convert_to_float(self):
        self.assertAlmostEqual(float(lit.RUPTURE_FORCE), 16.1, places=9)
        self.assertAlmostEqual(float(lit.CELL_LENGTH), 1.80, places=9)

    def test_formatting_includes_uncertainty_and_unit(self):
        self.assertEqual(lit.RUPTURE_FORCE.formatted, "16.1 ± 0.4 pN")

    def test_all_measurements_come_from_the_calibration_paper(self):
        """環境定数を除き、すべての値が同一の実測論文に由来すること。"""
        measured = [v for v in lit.ALL_VALUES if v.source != "tirado1979"]
        self.assertTrue(all(v.source == "costan2024" for v in measured))
        self.assertGreaterEqual(len(lit.values_from("costan2024")), 10)

    def test_parameters_table_lists_every_value(self):
        table = lit.parameters_markdown()
        for value in lit.ALL_VALUES:
            with self.subTest(symbol=value.symbol):
                self.assertIn(value.symbol, table)

    def test_bibliography_lists_every_reference(self):
        bibliography = lit.bibliography_markdown()
        for reference in lit.REFERENCES.values():
            self.assertIn(reference.doi, bibliography)


class AssumptionTests(unittest.TestCase):
    def test_ids_are_unique(self):
        ids = [item.id for item in ALL_ASSUMPTIONS]
        self.assertEqual(len(ids), len(set(ids)))

    def test_every_assumption_states_rationale_and_impact(self):
        for item in ALL_ASSUMPTIONS:
            with self.subTest(assumption=item.id):
                self.assertTrue(item.statement)
                self.assertTrue(item.rationale)
                self.assertTrue(item.impact)
                self.assertIsInstance(item.status, Status)

    def test_lookup_by_id_covers_every_assumption(self):
        self.assertEqual(len(BY_ID), len(ALL_ASSUMPTIONS))

    def test_uses_rejects_unknown_ids(self):
        with self.assertRaises(KeyError):
            uses("A1", "A999")

    def test_uses_returns_the_referenced_assumptions(self):
        selected = uses("A2", "A10")
        self.assertEqual([item.id for item in selected], ["A2", "A10"])

    def test_untested_assumptions_are_declared_as_such(self):
        """未検証の仮定が SUPPORTED に紛れていないこと。"""
        untested = {item.id for item in ALL_ASSUMPTIONS if item.status is Status.UNTESTED}
        self.assertIn("A6", untested)

    def test_documentation_renders_every_assumption(self):
        from adhesion.assumptions import assumptions_markdown

        document = assumptions_markdown()
        for item in ALL_ASSUMPTIONS:
            with self.subTest(assumption=item.id):
                self.assertIn(item.id, document)
                self.assertIn(item.rationale, document)


if __name__ == "__main__":
    unittest.main()
