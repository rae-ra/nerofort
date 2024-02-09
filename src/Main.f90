
program Main

    character (len=1) :: pathsep
    pathsep = "?"
#ifdef WOS
    pathsep = "\"
#endif
#ifdef UOS
    pathsep = "/"
#elif __unix
    pathsep = "/"
#endif
#if __APPLE__ || __MACH__
    pathsep = "/"
#endif

    Print*, "Hey"
    write (*, '( "pathsep is >", A1, "<")' )  pathsep



end program


