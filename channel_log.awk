BEGIN           { running= 0 }
/Test started/  { running= 1; next }
/Test finished/ { exit }
                { if(running) print }
