# Use mdbtools to export all tables in csv
# sudo apt-get install mdbtools

# Read schema
mdb-schema CDC_Obj2.accdb | 
# Find table names within that schema
grep -oP "CREATE TABLE \[\K.+(?=\])" | 
# Export each table into its own file under NIGGCWserol folder
while read line ; do mdb-export CDC_Obj2.accdb "$line" > "$line.csv" ; done
