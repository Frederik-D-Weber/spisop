@REM @CALL :sub >spisop.bat.log.txt
@REM @EXIT /b
@REM :sub
@COLOR 2F
@REM  Give the complete path to the toolbox folder (depends on your operating
@REM  system (Mac/Unix/Linux or Windows)

SET pathPrefix=D:\spisop_toolbox_beta2.3

@IF NOT exist %pathPrefix%\utility\convert_data\OpenBCIv3\convert_SD_file\SD_to_ascii_then_to_EEG\obci_convert_gain24.exe (
	echo spisop toolbox does not appear to be installed corecctly, please check if pathPrefix is set correctly!
	@EXIT /b
) 


SET settingsFile=split_pure_lindev_filter_home_jingyi_setup1_spind_slow_alpha_extension.txt
SET linearDeviationFile=linear_deviation_sleep_scoring.txt

@REM ECHO  convert edf_autoscale
@REM obci_convert_gain24.exe edf_autoscale %pathPrefix% %settingsFile% linearDeviationMontagePath=%linearDeviationFile%

ECHO  convert edf_0p1uVacc_cutoff
obci_convert_gain24.exe edf_0p1uVacc_cutoff %pathPrefix% %settingsFile% linearDeviationMontagePath=%linearDeviationFile%

@REM ECHO  convert edf_0p01uVacc_cutoff
@REM obci_convert_gain24.exe edf_0p01uVacc_cutoff %pathPrefix% %settingsFile% linearDeviationMontagePath=%linearDeviationFile%

ECHO  convert brainvision_int16
obci_convert_gain24.exe brainvision_int16 %pathPrefix% %settingsFile% linearDeviationMontagePath=%linearDeviationFile%

@REM ECHO  convert brainvision_int32
@REM obci_convert_gain24.exe brainvision_int32 %pathPrefix% %settingsFile% linearDeviationMontagePath=%linearDeviationFile%

@REM ECHO  convert brainvision_float32
@REM obci_convert_gain24.exe brainvision_float32 %pathPrefix% %settingsFile% linearDeviationMontagePath=%linearDeviationFile%

@COLOR

