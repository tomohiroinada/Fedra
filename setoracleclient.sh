#!/bin/sh

export FEDRA_PATH=/home/tomohiro/sim/FASER/fedra
	
### required for occi (access from c++):
export OCCIHOME=$FEDRA_PATH/instantclient
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCIHOME
export PATH=$PATH:$OCCIHOME

### required for $TNS_ADMIN/tnsnames.ora
export TNS_ADMIN=$FEDRA_PATH/oracle

### required for $ORACLE_HOME/sqplus/admin/glogon.sql  (optional)
export ORACLE_HOME=$FEDRA_PATH/oracle
