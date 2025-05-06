#!/bin/sh
doxygen Doxyfile
chcon -R -t httpd_sys_content_t html/
