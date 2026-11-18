"""The helper tools LTRquest fetches into --tools-dir, and what to say when it cannot.

detect, kmer2ltr and fastltr each clone a helper there on first use, and the
end-of-run clean-up removes a tools dir only when it holds nothing but those
helpers. Both facts live here because detect imports kmer2ltr and fastltr, so
neither of those can import them from detect.
"""

# What the fetchers put at the top of a tools dir: detect's TEsorter2, TRF-mod and
# sdust (and minimap2 and miniprot, when they are not on PATH), kmer2ltr's Kmer2LTR
# and fastltr's FASTAN, which FastLTR is built inside.
HELPERS = frozenset({"TEsorter2", "TRF-mod", "sdust", "minimap2", "miniprot",
                     "Kmer2LTR", "FASTAN"})

# Fetching a helper is the one step that needs the internet, and compute nodes often
# have none.
OFFLINE_HINT = ("Offline? Run once with internet using --tools-dir DIR, then pass the same "
                "--tools-dir; or use the container.")
