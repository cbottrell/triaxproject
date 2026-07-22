from pathlib import Path
import re
import unittest


README = Path(__file__).resolve().parents[1] / "README.md"


def _blockquote_content(line: str) -> str:
    """Return Markdown content after one blockquote marker, if present."""

    stripped = line.strip()
    if stripped.startswith(">"):
        return stripped[1:].lstrip()
    return stripped


class ReadmeMathTests(unittest.TestCase):
    def test_display_math_contains_no_setext_heading_markers(self):
        """A bare '=' or '-' line makes GitHub parse TeX as a heading."""

        in_display_math = False
        for line_number, line in enumerate(
            README.read_text(encoding="utf-8").splitlines(), start=1
        ):
            content = _blockquote_content(line)
            if content == "$$":
                in_display_math = not in_display_math
                continue

            if in_display_math:
                self.assertIsNone(
                    re.fullmatch(r"[=-]+", content),
                    f"Markdown heading marker inside display math on line {line_number}",
                )

        self.assertFalse(in_display_math, "README.md has an unclosed display-math block")

    def test_fragile_github_math_macros_are_absent(self):
        source = README.read_text(encoding="utf-8")
        forbidden = (r"\boxed", r"\!", r"\,", r"\mathbin", r"\operatorname")

        for macro in forbidden:
            with self.subTest(macro=macro):
                self.assertNotIn(macro, source)


if __name__ == "__main__":
    unittest.main()
