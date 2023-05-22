#! /usr/bin/env python3
#
# SNAP: Servere Nuclear Accident Programme
# Copyright (C) 1992-2017   Norwegian Meteorological Institute
# 
# This file is part of SNAP. SNAP is free software: you can 
# redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import sys
from PyQt5 import QtGui, QtWidgets
from Snappy.SnapController import SnapController
from Snappy.SnapControllerInverse import SnapControllerInverse
from Snappy.BrowserWidget import BrowserWidget
import Snappy.EEMEP.Resources
import Snappy.EEMEP.Controller

def main():
    app = QtWidgets.QApplication(sys.argv)
    tabs = QtWidgets.QTabWidget()
    snap = SnapController()
    snapBack = SnapControllerInverse()
    eemep = Snappy.EEMEP.Controller.Controller()
    tabs.addTab(snap.main, 'SNAP Nuclear Accident')
    tabs.addTab(snapBack.main, 'SNAP Backtracking')
    tabs.addTab(eemep.main, 'EEMEP Ash Cloud')
    tabs.resize(960,1024)
    tabs.setWindowTitle('SNAPpy')
    tabs.setWindowIcon(QtGui.QIcon(snap.res.getIconPath()))
    tabs.show()


    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
