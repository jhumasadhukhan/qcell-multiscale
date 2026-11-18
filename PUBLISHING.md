# Putting this on GitHub and getting a DOI

Two steps. The second is the one journals care about: a GitHub URL can change
or disappear, so Nature-family journals ask for a persistent identifier, and
the Code availability statement in the manuscript has a placeholder waiting
for it.

## 1. Create the repository

On github.com, create a new **public** repository named `qcell-multiscale`
under your account. Do not tick "add a README", "add .gitignore" or
"choose a license" — this folder already has all three, and the extra files
would collide on the first push.

Then, from inside this folder:

```bash
git init
git add .
git commit -m "QCell: quantum-to-cellular multiscale biosensor model"
git branch -M main
git remote add origin https://github.com/jhumasadhukhan/qcell-multiscale.git
git push -u origin main
```

If `git push` asks for a password, GitHub no longer accepts account passwords
over HTTPS. Either install the GitHub CLI and run `gh auth login`, or create a
personal access token (Settings → Developer settings → Personal access tokens)
and paste that in place of the password.

Check afterwards that `results/results.json` made it: it is about 320 kB and
is what lets a reader regenerate the figures without the 45-minute run. The
`.gitignore` excludes generated PNGs but deliberately keeps that file.

## 2. Mint a DOI with Zenodo

1. Sign in at [zenodo.org](https://zenodo.org) with your GitHub account.
2. Go to your Zenodo profile → **GitHub**, find `qcell-multiscale` in the list
   and switch the toggle **On**. This must happen *before* the release.
3. Back on GitHub: **Releases** → **Create a new release**, tag `v1.0.0`,
   title "Initial release accompanying the manuscript", publish it.
4. Zenodo picks the release up within a few minutes and issues two DOIs: one
   for `v1.0.0` and a **concept DOI** that always resolves to the newest
   version. **Cite the concept DOI** in the paper.
5. Put that DOI in two places:
   - the manuscript's Code availability section, replacing
     `[Zenodo DOI to be inserted on acceptance]`
   - the `identifiers:` field of `CITATION.cff` in this repository

Zenodo reads `CITATION.cff`, so the author name, affiliation and title in the
archive record come from that file. Check it says what you want before you cut
the release.

## A note on timing

Some journals want the repository public at submission; others accept a
private repository plus a reviewer link until acceptance. Nature Chemical
Engineering expects code to be available to referees at submission, so public
from the start is the simpler path. If you would rather keep it private for
now, create it private and add a note in the cover letter saying the DOI will
be minted on acceptance — but remember that the Zenodo toggle only sees public
repositories.
