# put col1 to end of line

{for (i = 2; i<= NF; i++)
  printf("%s\t",$i);
printf("%s\n",$1);
}
