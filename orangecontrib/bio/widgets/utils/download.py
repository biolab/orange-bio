import itertools

from PyQt4.QtCore import pyqtSignal as Signal
from Orange.utils import serverfiles

from Orange.OrangeWidgets.OWConcurrent import Task


class EnsureDownloaded(Task):
    advance = Signal()
    progress = Signal(float)

    def __init__(self, files_list, parent=None):
        Task.__init__(self, parent)
        self.files_list = files_list

    def run(self):
        nfiles = len(self.files_list)
        count = 100 * nfiles
        counter = itertools.count()

        def advance():
            self.advance.emit()
            self.progress.emit(100.0 * next(counter) / count)

        for domain, filename in self.files_list:
            ensure_downloaded(domain, filename, advance=advance)


def ensure_downloaded(domain, filename, advance=None):
    serverfiles.localpath_download(domain, filename, callback=advance)
