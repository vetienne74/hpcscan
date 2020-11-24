#!/usr/bin/env bash

root=$1
if test -z "$root" ; then
    root=..
fi

git log $root | head -n3 > ./version.work
commit=$(grep -i commit ./version.work)
author=$(grep -i author ./version.work)
date=$(grep -i date ./version.work)
date_compile=$(date +%s)

rm ./version.work

echo "const char HPCSCAN_GIT_COMMIT[] = \"$commit\" ;"
echo "const char HPCSCAN_GIT_AUTHOR[] = \"$author\" ;"
echo "const char HPCSCAN_GIT_DATE[] = \"$date\" ;"
echo "const time_t HPCSCAN_COMPILE_DATE = $date_compile ;"
