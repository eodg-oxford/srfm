# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
import os
import re
import warnings

project = 'SRFM'
copyright = '2026, Antonin Knizek, Roy Grainger'
author = 'Antonin Knizek, Roy Grainger'
release = '0.0.2'

sys.path.insert(0, os.path.abspath('../../src'))


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.napoleon",
              "sphinx.ext.autodoc",
              "sphinx.ext.apidoc",
              "sphinx.ext.todo"
              ]

templates_path = ['_templates']
exclude_patterns = []
napoleon_custom_sections = [('Returns', 'params_style')]

warnings.filterwarnings("ignore", category=SyntaxWarning)

_CODE_FENCE_RE = re.compile(r"```(\w+)?\s*\n(.*?)```", re.DOTALL)


def _format_code_block(match: re.Match) -> str:
    """Convert Markdown-style code fences to reStructuredText blocks."""
    language = (match.group(1) or "text").strip() or "text"
    body = match.group(2).strip("\n")
    indented = "\n".join(f"    {line}" if line else "    " for line in body.splitlines())
    return f"\n.. code-block:: {language}\n\n{indented}\n"


def _fix_definition_list_spacing(lines: list[str]) -> list[str]:
    """Ensure definition lists are separated from following paragraphs."""
    fixed: list[str] = []
    prev_was_definition = False
    for line in lines:
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())
        is_blank = stripped == ""
        is_definition = indent == 2 and ":" in stripped
        is_paragraph = indent == 2 and not is_definition and not is_blank
        if is_paragraph and prev_was_definition and (not fixed or fixed[-1].strip()):
            fixed.append("")
            prev_was_definition = False
        fixed.append(line)
        if is_definition:
            prev_was_definition = True
        elif not is_blank and indent <= 2:
            prev_was_definition = False
    if prev_was_definition and (not fixed or fixed[-1].strip()):
        fixed.append("")
    return fixed


def _collapse_attribute_note_block(lines: list[str]) -> list[str]:
    """Replace malformed attribute directives with a short note."""
    target = ".. attribute:: A range of quality/flag fields"
    result: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(target):
            collected: list[str] = []
            while i < len(lines) and lines[i].startswith(".. attribute::"):
                collected.append(lines[i].split("::", 1)[1].strip())
                i += 1
                while i < len(lines) and (lines[i].startswith("   ") or not lines[i].strip()):
                    i += 1
            merged = " ".join(part for part in collected if part)
            note_body = merged or (
                "Additional quality flag fields are initialised in `_inivar`."
            )
            result.extend([".. note::", "", f"   {note_body}", ""])
            continue
        result.append(line)
        i += 1
    return result


def _ensure_param_type_spacing(lines: list[str]) -> list[str]:
    """Insert blank lines between parameter/type field-list entries."""
    fixed: list[str] = []
    seen_param = False
    for line in lines:
        if line.startswith(":param"):
            if seen_param and (not fixed or fixed[-1].strip()):
                fixed.append("")
            seen_param = True
        fixed.append(line)
        if not line.strip():
            seen_param = False
    return fixed


def _ensure_bullet_spacing(lines: list[str]) -> list[str]:
    """Ensure lists in narrative sections are separated from following text."""
    fixed: list[str] = []
    prev_was_bullet = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Reference:") and prev_was_bullet and (not fixed or fixed[-1].strip()):
            fixed.append("")
        fixed.append(line)
        if stripped.startswith("* "):
            prev_was_bullet = True
        elif not stripped:
            prev_was_bullet = False
    return fixed


def _cleanup_docstring(_app, _what, _name, _obj, _options, lines):
    """Rewrite problematic patterns before Sphinx parses the docstring."""
    raw = "\n".join(lines)
    cleaned = _CODE_FENCE_RE.sub(_format_code_block, raw)
    working = cleaned.splitlines()
    fixed = _fix_definition_list_spacing(working)
    rewritten = _collapse_attribute_note_block(fixed)
    spaced = _ensure_param_type_spacing(rewritten)
    separated = _ensure_bullet_spacing(spaced)
    lines[:] = separated



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']


def setup(app):
    app.connect("autodoc-process-docstring", _cleanup_docstring, priority=1000)
