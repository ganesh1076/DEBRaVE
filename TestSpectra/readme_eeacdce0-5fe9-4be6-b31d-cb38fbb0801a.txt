This file contains an overview of the structure and content of your download request from the ESO Science Archive.


For every downloaded dataset, files are listed below with the following structure:

dataset_name
        - archive_file_name (technical name, as saved on disk)	original_file_name (user name, contains relevant information) category size


Please note that, depending on your operating system and method of download, at download time the colons (:) in the archive_file_name as listed below may be replaced by underscores (_).


In order to rename the files on disk from the technical archive_file_name to the more meaningful original_file_name, run the following shell command:
    cat THIS_FILE | awk '$2 ~ /^ADP/ {print "test -f",$2,"&& mv",$2,$3}' | sh


In case you have requested cutouts, the file name on disk contains the TARGET name that you have provided as input. To order files by it when listing them, run the following shell command:
    cat THIS_FILE | awk '$2 ~ /^ADP/ {print $2}' | sort -t_ -k3,3


Your feedback regarding the data quality of the downloaded data products is greatly appreciated. Please contact the ESO Archive Science Group via https://support.eso.org/ , subject: Phase 3 ... thanks!

The download includes contributions from the following collections:
Ref(0)	HARPS	https://doi.eso.org/10.18727/archive/33	harps_20230215.pdf	https://www.eso.org/rm/api/v1/public/releaseDescriptions/72

Publications based on observations collected at ESO telescopes must acknowledge this fact (please see: http://archive.eso.org/cms/eso-data-access-policy.html#acknowledgement). In particular, please include a reference to the corresponding DOI(s). They are listed in the third column in the table above and referenced below for each dataset. The following shell command lists them:

	cat THIS_FILE | awk -F/ '$1 ~ /^Ref\(/ {print $0,$NF}' | awk '{print $2, $3}' | sort | uniq


Each collection is described in detail in the corresponding Release Description. They can be downloaded with the following shell command:

	cat THIS_FILE | awk -F/ '$1 ~ /^Ref\(/ {print $0,$NF}' | awk '{printf("curl -o %s_%s %s\n", $6, $4, $5)}' | sh

ADP.2014-09-26T16:53:10.623 Ref(0)
	- ADP.2014-09-26T16:53:10.623.fits	HARPS.2012-09-08T09:47:59.185_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-26T16:52:22.940 Ref(0)
	- ADP.2014-09-26T16:52:22.940.fits	HARPS.2012-09-08T08:34:21.289_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-25T15:37:21.833 Ref(0)
	- ADP.2014-09-25T15:37:21.833.fits	HARPS.2011-09-08T09:24:29.528_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-25T15:35:09.290 Ref(0)
	- ADP.2014-09-25T15:35:09.290.fits	HARPS.2011-09-07T08:53:03.514_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-25T15:33:52.153 Ref(0)
	- ADP.2014-09-25T15:33:52.153.fits	HARPS.2011-09-09T08:06:33.362_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-25T15:35:55.137 Ref(0)
	- ADP.2014-09-25T15:35:55.137.fits	HARPS.2011-09-07T07:45:25.731_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-26T16:51:28.527 Ref(0)
	- ADP.2014-09-26T16:51:28.527.fits	HARPS.2012-07-30T09:46:48.610_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-26T16:51:08.790 Ref(0)
	- ADP.2014-09-26T16:51:08.790.fits	HARPS.2012-07-29T10:23:27.539_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-26T16:52:47.840 Ref(0)
	- ADP.2014-09-26T16:52:47.840.fits	HARPS.2012-09-09T09:48:25.905_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
ADP.2014-09-25T15:36:23.730 Ref(0)
	- ADP.2014-09-25T15:36:23.730.fits	HARPS.2011-09-08T07:35:37.261_s1d_A_DRS_HARPS_3.5_ESOSDP.fits	SCIENCE.SPECTRUM	5264640
