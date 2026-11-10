"""Small, tool-independent helpers for human and machine-readable reports."""
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def jsonable(value):
    if isinstance(value, Path):
        return str(value.expanduser().absolute())
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def atomic_text(path, text):
    path = Path(path)
    tmp = path.with_name(f'.{path.name}.{uuid4().hex}.tmp')
    try:
        tmp.write_text(text, encoding='utf-8')
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)
