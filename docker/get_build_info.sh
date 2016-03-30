#!/usr/bin/env bash

if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then
    echo "TODO"
else
    export PULL_REQUEST_TITLE=$(curl -s "https://api.github.com/repos/ORNL-CEES/Cap/pulls/${TRAVIS_PULL_REQUEST}" | jq '.title' | sed 's/"//g')
    echo "Pull Request #${TRAVIS_PULL_REQUEST}: ${PULL_REQUEST_TITLE}"
fi

