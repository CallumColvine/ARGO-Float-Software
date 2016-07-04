rem call %cd%\ui_Files\build.bat
rem call %cd%\ui_Classes\run.bat
@ECHO off

ECHO Rebuild ui files
pyside-uic %cd%\ui_Files\MainApp.ui -o %cd%\ui_Files\ui_MainApp.py
pyside-uic %cd%\ui_Files\timeseriesapp.ui -o %cd%\ui_Files\ui_timeseriesapp.py
pyside-uic %cd%\ui_Files\circulationapp.ui -o %cd%\ui_Files\ui_circulationapp.py
ECHO ui File build complete

ECHO Run the GUI
REM pyside-uic .\ui_Files\MainApp.ui -o ui_MainApp.py
ECHO ----- End of Batch File. Running GUI now -----
python %cd%\ui_Classes\MainApp.py 
