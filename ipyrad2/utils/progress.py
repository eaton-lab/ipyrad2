import datetime
import sys
import time


class ProgressBar(object):
    """
    Print pretty progress bar
    """
    def __init__(self, njobs, start=None, message=""):
        self.njobs = njobs
        self.start = (start if start else time.time())
        self.message = message
        self.finished = 0

    @property
    def progress(self):
        return 100 * (self.finished / float(self.njobs))

    @property
    def elapsed(self):
        return datetime.timedelta(seconds=int(time.time() - self.start))

    def update(self):
        # build the bar
        hashes = '#' * int(self.progress / 5.)
        nohash = ' ' * int(20 - len(hashes))

        # print to stderr
        print("\r[{}] {:>3}% {} | {:<12} ".format(*[
            hashes + nohash,
            int(self.progress),
            self.elapsed,
            self.message,
        ]), end="")
        sys.stdout.flush()
