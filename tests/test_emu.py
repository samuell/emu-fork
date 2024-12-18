import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import emu


def test_generate_outfile_path():
    test_cases = [
        {
            "output_dir": "results",
            "output_basename": "",
            "input_file": "some_file.fastq",
            "want_outfile": "results/some_file.fastq",
        },
        {
            "output_dir": "results",
            "output_basename": "",
            "input_file": ["some_file.fastq", "another_file.fastq"],
            "want_outfile": "results/some_file-another_file",
        },
        {
            "output_dir": "results",
            "output_basename": "my-output-basename",
            "input_file": "some_file.fastq",
            "want_outfile": "results/my-output-basename",
        },
    ]

    for tc in test_cases:
        have_outfile = emu.create_outfile_path(
            tc["output_dir"], tc["output_basename"], tc["input_file"]
        )
        assert have_outfile == tc["want_outfile"]
