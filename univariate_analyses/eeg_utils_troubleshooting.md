# Installing gfortran on macOS for R package `eegUtils`

If the installation is conflicting with software from FSL (eg, you see
"fsl" in the filepath in error messages), run the following in your current
R session:


``` R
Sys.setenv(PATH = gsub("(/usr/local/fsl/bin:|:/usr/local/fsl/bin)", "", Sys.getenv("PATH")))
Sys.unsetenv("CPLUS_INCLUDE_PATH")
Sys.unsetenv("CPPFLAGS")
```

and try installing again:

```R
remotes::install_github("craddm/eegUtils")
```



`eegUtils` depends on compiled Fortran code, and macOS users sometimes hit
issues
when the installer can't find `gfortran` in `/opt/gfortran`.


If installing with:

``` r
remotes::install_github('craddm/eegUtils')
```

fails due to a missing gfortran version, follow these steps:

**Check for and uninstall brew gfortran**

Type `which gfortran` in your terminal. If you see `brew` in the file path,
the wrong version is installed. Run:

```bash
brew remove gfortran
```

**Install correct version**

1. Visit: https://github.com/fxcoudert/gfortran-for-macOS/releases
2. Find the `.dmg` file for the version it's looking for (eg,
   `12.2-arm64.dmg` for Apple Silicon)
3. Right-click and **Copy Link Address** for that `.dmg`
4. In your terminal, run (replacing the URL below with your copied link):

```bash
curl -L -o gfortran.dmg https://github.com/fxcoudert/gfortran-for-macOS/releases/download/v12.2/gfortran-12.2-arm64.dmg
hdiutil attach gfortran.dmg
cd /Volumes/gfortran-* # Just cd into whatever version gets installed here
sudo installer -pkg *.pkg -target /
cd ~
hdiutil detach /Volumes/gfortran-*
```

5. **Move gfortran to correct location.** The gfortran library needs to be
   in the `/opt` directory in order for eegUtils to find it. Find where
   gfortran got installed to with 

```bash
ls -l | which gfortran
```
Then move it to the correct location (eg, `sudo mv /usr/local/gfortran
/opt/gfortran`).

6. Try installing `eegUtils` again:

```r
remotes::install_github('craddm/eegUtils')
```

7. If it still crashes and mentions it "can't find directory
   `/opt/gfortran/lib/gcc/aarch64-apple-darwinXX.YY`, manually rename the
   folder in `/opt/gfortran/lib/gcc/` to match the version it's requesting.


