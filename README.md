# rt_lookup
look up functions for adding retention time in msp files

the default msp dir on hive is: /quobyte/metabolomicsgrp/fanzhou/raw_msp

the default prediction dir on hive is:/quobyte/metabolomicsgrp/fanzhou/predictions (for look up tables)


The msp files with predicted RT added will also be saved at default msp dir.

sample usage of the script:

python add_rt_msp.py --library_name nist23 --column_name rp --msp_dir /quobyte/metabolomicsgrp/fanzhou/raw_msp --prediction_dir /quobyte/metabolomicsgrp/fanzhou/predictions

if used on hive, the --msp_dir and --prediction_dir can be omitted since they are set as default

