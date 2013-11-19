# channel_log.awk
#
# Extract samples from test results
#
# This code is experimental, and error-handling is primitive.
#
# @LICENSE(NICTA)

BEGIN           { running= 0 }
/Test started/  { running= 1; next }
/Test finished/ { exit }
                { if(running) print }
