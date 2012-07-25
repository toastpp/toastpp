import os
import sys

major = "%d" % sys.version_info[0]
minor = "%d" % sys.version_info[1]
pyver = "python" + major + "." + minor
pytoastdir = os.path.expandvars("$TOASTVER/lib/" + pyver + "/site-packages")
sys.path.append(pytoastdir)
