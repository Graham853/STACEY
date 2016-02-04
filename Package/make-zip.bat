SET ZIPFILE=STACEY.addon.v1.1.1.zip
SET SEVENZIP="C:\Program Files\7-Zip\7z.exe"

DEL /P %ZIPFILE%

%SEVENZIP% a %ZIPFILE% doc\*
%SEVENZIP% a %ZIPFILE% examples\*
%SEVENZIP% a %ZIPFILE% lib\*
%SEVENZIP% a %ZIPFILE% templates\*
%SEVENZIP% a %ZIPFILE% stacey.src.jar
%SEVENZIP% a %ZIPFILE% version.xml
pause


