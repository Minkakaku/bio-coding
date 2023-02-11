# Function get from https://stackoverflow.com/questions/10674992/convert-a-character-vector-of-mixed-numbers-fractions-and-integers-to-numeric?rq=1
# With little modifications
mixedToFloat <- function(x) {
    x <- sapply(x, as.character)
    is.integer <- grepl("^-?\\d+$", x)
    is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
    is.float <- grepl("^-?\\d+\\.\\d+$", x)
    is.mixed <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
    stopifnot(all(is.integer | is.fraction | is.float | is.mixed))

    numbers <- strsplit(x, "[ /]")

    ifelse(is.integer, as.numeric(sapply(numbers, `[`, 1)),
        ifelse(is.float, as.numeric(sapply(numbers, `[`, 1)),
            ifelse(is.fraction, as.numeric(sapply(numbers, `[`, 1)) /
                as.numeric(sapply(numbers, `[`, 2)),
            as.numeric(sapply(numbers, `[`, 1)) +
                as.numeric(sapply(numbers, `[`, 2)) /
                    as.numeric(sapply(numbers, `[`, 3))
            )
        )
    )
}