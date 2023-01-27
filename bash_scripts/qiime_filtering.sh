#!/bin/bash

# filter sequences based on the sample

qiime feature-table filter-features \
  --i-table table.qza
  --o-filtered-table <output_name>
