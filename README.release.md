# Release guidelines

This recipe might need improvements, but keep the spirit in mind and update as
appropriate.

* Make sure all documentation and references are updated and all tests are
  passed.
* Update `Release Notes.md`.
* Bump version (in setup.py and doc/conf.py) to reflect the new release
* Copy current `develop` to `release`:
```
git checkout develop@{0}  # get working tree from "development", detach HEAD
git reset --soft stable   # reposition detached HEAD on "stable"
git commit                # enter the appropriate commit message
git branch temp           # create a temporary branch "temp" at HEAD
git checkout temp         # get on the new temporary branch
git branch -M release     # rename "temp" to "stable"
```
* Tag the final release as vx.y or vx.y.z or, possibly, create a vx.y branch for
  the release.
* go back to develop `git checkout develop`
* Set the version (in setup.py) to the next development-cycle version number,
  e.g., `x.u-dev`.


## Publish release version

Push the release branch(es) and tag to github
```
git push github release
git push github vx.y
```


## Create github release

Go the the github project page and create a new release from the pushed tag
vx.y.

Provide a release title "Public release vx.y" and add the release notes from
`Release Notes.md` as description of the github release.


<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 80
End:
-->
