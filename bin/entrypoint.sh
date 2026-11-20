#!/bin/sh
# Entry point for the LTR_peel container image.
#
# Apptainer turns a Docker ENTRYPOINT into the SIF's runscript, which is what
# decides whether `./ltr_peel.sif --genome x.fa` means anything. An image with
# only a CMD gets a runscript that *replaces* the command with whatever the user
# typed, so 1.0.1 answered
#
#     ./ltr_peel.sif --help
#     FATAL: "--help": executable file not found in $PATH
#
# Setting an ENTRYPOINT makes Apptainer prepend it to the user's arguments
# instead, and this script decides what those arguments mean:
#
#     ./ltr_peel.sif --genome x.fa       -> ltr_peel --genome x.fa
#     ./ltr_peel.sif ltr_peel-gff3 -h    -> ltr_peel-gff3 -h
#     docker run IMG ltr_peel --help     -> ltr_peel --help      (as documented)
#     docker run IMG /bin/bash -c '...'  -> /bin/bash -c '...'   (Nextflow)
#     apptainer exec IMG ltr_peel ...    -> never reaches here; exec skips it
#
# The rule: a first argument that names a runnable command is a command;
# anything else is an ltr_peel argument. That is strictly more permissive than
# the empty ENTRYPOINT it replaces, so no invocation that worked before stops
# working.

set -eu

case "${1-}" in
    ""|-*)
        # No arguments, or the first one is a flag: both belong to ltr_peel.
        # Matched before `command -v` because `command -v --help` is not
        # portable across shells.
        ;;
    *)
        if command -v "$1" > /dev/null 2>&1; then
            exec "$@"
        fi
        ;;
esac

exec ltr_peel "$@"
