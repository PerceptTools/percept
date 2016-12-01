#!/bin/bash
./mesh_adapt --version
if (($? > 0))
then
  echo "built in Sierra"
else
  echo "not built in Sierra"
fi
