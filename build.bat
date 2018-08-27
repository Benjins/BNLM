@echo off

cl /Od /Zi /DBNS_DEBUG test_main.cpp /Fetest_main.exe

::g++ -std=c++11 -O0 -g -DBNS_DEBUG test_main.cpp -o test_main.exe
