ECHO Rebuild ui files
@ECHO off
pyside-uic %cd%\MainApp.ui -o %cd%\ui_MainApp.py
pyside-uic %cd%\timeseriesapp.ui -o %cd%\ui_timeseriesapp.py
pyside-uic %cd%\circulationapp.ui -o %cd%\ui_circulationapp.py
ECHO ui File build complete
