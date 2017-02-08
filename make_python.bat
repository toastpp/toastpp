SET VS90COMNTOOLS=%VS140COMNTOOLS%
"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\x86_amd64\vcvarsx86_amd64.bat"
pip install numpy
cd src\python
python setup.py build --build-base=../../win/x64/Release/python
python setup.py install
cd ..\..
