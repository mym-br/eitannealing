from typing import TextIO


class ValidationError(Exception):
    pass


class InvalidTripletCountError(ValidationError):
    def __init__(self, element_count):
        super().__init__(
            f"The number of elements in the file ({element_count}) is not a multiple of 3."
        )
        self.element_count = element_count


class InvalidElectrodeIndexError(ValidationError):
    def __init__(self, index, electrodes_count):
        super().__init__(
            f"Electrode index {index} is outside the range [1, {electrodes_count}]."
        )
        self.index = index
        self.electrodes_count = electrodes_count


class InvalidElectrodesCountError(ValidationError):
    def __init__(self, injections_count, electrodes_count):
        super().__init__(
            f"The number of injections in the file ({injections_count}) is not equal to the electrodes count ({electrodes_count})."
        )
        self.injections_count = injections_count
        self.electrodes_count = electrodes_count


class InvalidFloatValueError(ValidationError):
    def __init__(self, index, value):
        super().__init__(f"Could not parse element at index {index}: {value}")
        self.index = index
        self.value = value


def validate_electrode_index(index: int, electrodes_count: int) -> None:
    if index < 1 or index > electrodes_count:
        raise InvalidElectrodeIndexError(index, electrodes_count)


def validate_currents(currents_file_content: TextIO, electrodes_count: int) -> None:
    content = currents_file_content.read().strip()

    # Split the content into individual elements
    elements = content.split()

    if len(elements) % 3 != 0:
        raise InvalidTripletCountError(len(elements))

    if len(elements) / 3 != electrodes_count:
        raise InvalidElectrodesCountError(len(elements) // 3, electrodes_count)

    for i in range(0, len(elements), 3):
        idx = i
        try:
            validate_electrode_index(int(elements[idx]), electrodes_count)
            validate_electrode_index(int(elements[idx := i + 1]), electrodes_count)
            float(elements[idx := i + 2])
        except ValueError as exc:
            raise InvalidFloatValueError(idx, elements[idx]) from exc
