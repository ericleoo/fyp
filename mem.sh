free | grep Mem | awk '{printf("free: %.4f%\t%d\nused: %.4f%\t%d\n", $7/$2 * 100.0, $7/1024, $3/$2 * 100.0, $3/1024)}'
#free | grep Mem | awk '{printf("%d\n",$7)}'
