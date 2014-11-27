For /D %%i in (Job*) do rmdir "%%i" /S /Q
For %%i in (Job*.mat) do del "%%i" /S /Q
del matlab_metadata.mat
del misc_Info.mat
