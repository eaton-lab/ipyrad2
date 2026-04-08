import sys
import time

from loguru import logger

from .logger import is_log_level_enabled


class ProgressBar(object):
    """
    Print pretty progress bar
    """
    def __init__(self, njobs, start=None, message=""):
        self.njobs = njobs
        self.start = (start if start else time.time())
        self.message = message
        self.finished = 0
        self._visible = False
        self._last_render = None

    @property
    def progress(self):
        return 100 * (self.finished / float(self.njobs))

    def render(self) -> str:
        """Return the current progress line without logger formatting."""
        hashes = '#' * int(self.progress / 5.)
        nohash = ' ' * int(20 - len(hashes))
        return "[{}] {:>3}% | {:<12} ".format(*[
            hashes + nohash,
            int(self.progress),
            self.message,
        ])

    def update(self):
        """Emit the current progress bar at INFO level with log-style prefix."""
        if not is_log_level_enabled("INFO"):
            return
        rendered = self.render()
        if rendered == self._last_render:
            return
        self._last_render = rendered
        self._visible = True
        logger.bind(end="\r").opt(depth=1).info(rendered)

    def close(self) -> None:
        """Finish the transient progress bar line if it was visible."""
        if not self._visible:
            return
        sys.stderr.write("\n")
        sys.stderr.flush()
        self._visible = False
        self._last_render = None
