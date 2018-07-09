#!/bin/sh

. ~/.bashrc

echo 'MakeGAPDocDoc("./","orders.xml",[],"orders");' | gap
