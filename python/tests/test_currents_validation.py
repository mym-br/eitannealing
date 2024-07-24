import os
import tempfile
import unittest

from pyeitsolver import EitComplexSolver, EitSolver
from pyeitsolver.validation import (
    InvalidElectrodeIndexError,
    InvalidElectrodesCountError,
    InvalidFloatValueError,
    InvalidTripletCountError,
    ValidationError,
)

# Example valid input
valid_input = "1 5 1.9020186995339855e-03 2 6 1.9020186995339855e-03 3 7 1.9020186995339855e-03 4 8 1.9020186995339855e-03 5 9 1.9020186995339855e-03 6 10 1.9020186995339855e-03 7 11 1.9020186995339855e-03 8 12 1.9020186995339855e-03 9 13 1.9020186995339855e-03 10 14 1.9020186995339855e-03 11 15 1.9020186995339855e-03 12 16 1.9020186995339855e-03 13 17 1.9020186995339855e-03 14 18 1.9020186995339855e-03 15 19 1.9020186995339855e-03 16 20 1.9020186995339855e-03 17 21 1.9020186995339855e-03 18 22 1.9020186995339855e-03 19 23 1.9020186995339855e-03 20 24 1.9020186995339855e-03 21 25 1.9020186995339855e-03 22 26 1.9020186995339855e-03 23 27 1.9020186995339855e-03 24 28 1.9020186995339855e-03 25 29 1.9020186995339855e-03 26 30 1.9020186995339855e-03 27 31 1.9020186995339855e-03 28 32 1.9020186995339855e-03 29 1 1.9020186995339855e-03 30 2 1.9020186995339855e-03 31 3 1.9020186995339855e-03 32 4 1.9020186995339855e-03 "


class TestValidateCurrents(unittest.TestCase):
    def setUp(self):
        self.electrodes_count = 32
        self.valid_input = valid_input
        self.mesh_filename = os.path.join(
            os.path.dirname(__file__), "../../data/circular_A_2D.msh"
        )

    def construct_solver(self, currents_file_content: str, currents_count: int = 32):
        with tempfile.NamedTemporaryFile(mode="w", delete_on_close=False) as fp:
            fp.write(currents_file_content)
            fp.close()
            EitSolver(self.mesh_filename, fp.name, currents_count)

    def test_valid_input(self):
        try:
            self.construct_solver(self.valid_input)
        except ValidationError:
            self.fail("validate_currents raised ValidationError unexpectedly!")

    def test_invalid_triplet_count(self):
        invalid_input = self.valid_input + "33"
        with self.assertRaises(InvalidTripletCountError) as context:
            self.construct_solver(invalid_input)
        self.assertEqual(context.exception.element_count, 97)

    def test_invalid_electrode_index(self):
        invalid_input = "0 5 " + self.valid_input[4:]
        with self.assertRaises(InvalidElectrodeIndexError) as context:
            self.construct_solver(invalid_input)
        self.assertEqual(context.exception.index, 0)
        self.assertEqual(context.exception.electrodes_count, 32)

        invalid_input = "1 33 " + self.valid_input[4:]
        with self.assertRaises(InvalidElectrodeIndexError) as context:
            self.construct_solver(invalid_input)
        self.assertEqual(context.exception.index, 33)
        self.assertEqual(context.exception.electrodes_count, 32)

    def test_invalid_float_value(self):
        invalid_input = "1 5 abc" + self.valid_input[26:]
        with self.assertRaises(InvalidFloatValueError) as context:
            self.construct_solver(invalid_input)
        self.assertEqual(context.exception.index, 2)
        self.assertEqual(context.exception.value, "abc")

    def test_invalid_electrodes_count(self):
        with self.assertRaises(InvalidElectrodesCountError) as context:
            self.construct_solver(self.valid_input, self.electrodes_count - 1)
        self.assertEqual(context.exception.injections_count, self.electrodes_count)
        self.assertEqual(context.exception.electrodes_count, self.electrodes_count - 1)


class TestValidateCurrentsComplex(TestValidateCurrents):

    def construct_solver(self, currents_file_content: str, currents_count: int = 32):
        with tempfile.NamedTemporaryFile(mode="w", delete_on_close=False) as fp:
            fp.write(currents_file_content)
            fp.close()
            EitComplexSolver(self.mesh_filename, fp.name, currents_count)


if __name__ == "__main__":
    unittest.main()
