.First.lib <- function(lib, pkg)
{
    library.dynam("gss", pkg, lib)
}
project <- function (object,...)
{
    UseMethod("project")
}
