#!/usr/bin/env bash

if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then
    curl -s "https://api.travis-ci.org/repos/${TRAVIS_REPO_SLUG}/builds/${TRAVIS_BUILD_ID}" | jq '.message' | sed -e 's/\n.*//' -e 's/"//g'
else
    export PULL_REQUEST_TITLE=$(curl -s "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/pulls/${TRAVIS_PULL_REQUEST}" | jq '.title' | sed 's/"//g')
    echo "PR #${TRAVIS_PULL_REQUEST}: ${PULL_REQUEST_TITLE}"
fi

