#!/bin/bash
for d in "$@"; do
    rm -rf ${d}/*out*
    rm -rf ${d}/*time*
    rm -rf ${d}/sal_*
    rm -rf ${d}/*eval*
    rm -rf ${d}/mcr16*
    rm -rf ${d}/concordance.json
    rm -rf /tmp/jcascitt/si/${d}
done
