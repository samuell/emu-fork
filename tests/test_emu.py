import os
import pytest
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import emu


@pytest.mark.parametrize(
    "output_dir, output_basename, input_file, expected_output_basepath",
    [
        (
            "results",  # output_dir
            "",  # output_basename
            "some_file.fastq",  # input_file
            "results/some_file.fastq",  # expected_output_basepath
        ),
        (
            "results",  # similar as above ...
            "",
            ["some_file.fastq", "another_file.fastq"],
            "results/some_file-another_file",
        ),
        (
            "results",
            "my-output-basename",
            "some_file.fastq",
            "results/my-output-basename",
        ),
    ],
)
def test_generate_outfile_path(
    output_dir, output_basename, input_file, expected_output_basepath
):
    got_output_basepath = emu.create_output_basepath(
        output_dir, output_basename, input_file
    )
    assert got_output_basepath == expected_output_basepath
