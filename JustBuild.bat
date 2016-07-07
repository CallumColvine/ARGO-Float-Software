@ECHO off

ECHO Rebuild ui files
pyside-uic %cd%\ui_Files\MainApp.ui -o %cd%\ui_Files\ui_MainApp.py
pyside-uic %cd%\ui_Files\timeseriesapp.ui -o %cd%\ui_Files\ui_timeseriesapp.py
pyside-uic %cd%\ui_Files\circulationapp.ui -o %cd%\ui_Files\ui_circulationapp.py
ECHO ui File build complete
