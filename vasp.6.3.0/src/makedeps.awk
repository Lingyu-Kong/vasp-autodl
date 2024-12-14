FNR==1{
    file = FILENAME
    sub(/.(f|F|fpp|FPP|for|FOR|ftn|FTN|f90|F90|f95|F95|f03|F03|f08|F08)$/,"",file)
    files[file]=file
}
tolower($1) == "module" && tolower($0) !~ /^[^!]+(subroutine|function|procedure)[[:blank:]]+[^!]/{
    name = tolower($2)
    sub(/!.*$/,"",name)
    whereIsMod[name]=file
}
tolower($1) == "use"{
    name = tolower($0)
    sub(/^[ \t]*use[ \t]*/,"",name)
    sub(/^(.*::)?[ \t]*/,"",name)
    sub(/[ \t]*((,|!).*)?$/,"",name)
    numUses[file]++
    usedMods[numUses[file],file] = name
}
END{
    for (i in files){
        printf("%s.o : %s.f90",files[i],i)
        for (j=1;j<=numUses[i];j++){
            if(!alreadyListed[whereIsMod[usedMods[j,i]],files[i]] &&  whereIsMod[usedMods[j,i]] != ""){
                alreadyListed[whereIsMod[usedMods[j,i]],files[i]] = 1
                if (whereIsMod[usedMods[j,i]] != files[i]){
                    printf(" %s.o",whereIsMod[usedMods[j,i]])
                }
            }
        }
        printf("\n")
    }
    if (numUses["version"] != 1 || usedMods[1,"version"] != "string"){
        printf("Your awk does not seem to work as expected.\n")
        exit 1
    }
}
