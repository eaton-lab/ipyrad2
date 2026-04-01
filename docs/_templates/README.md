# Docs Authoring Guide

This directory holds reusable draft templates for the ipyrad2 documentation site.

## Which Template to Use

- `section-landing.md`: use for section hubs that orient the reader and link to child pages.
- `concept-guide.md`: use for narrative pages that explain ideas, workflows, or project philosophy.
- `command-reference.md`: use for command pages and tool references.

## Naming and Layout Conventions

- Keep the published docs tree under `docs/`.
- Use `index.md` for section landing pages.
- Use lowercase kebab-case for page filenames.
- Keep templates in `docs/_templates/` and out of navigation.

## Markdown-First Rule

- Write docs in Markdown by default.
- Keep examples executable where possible.
- Prefer concise placeholders over long unfinished prose in draft pages.

## Notebook Exception Policy

- Notebooks are a later exception, not the default authoring path.
- When a notebook is eventually added, it should produce or pair with a published Markdown page.
- Do not treat raw notebooks as the canonical published documentation tree.
