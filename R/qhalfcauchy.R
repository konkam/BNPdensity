qhalfcauchy <-
function (p, location = 0, scale = 1) 
{
    qcauchy(p * (1 - pcauchy(0, location, scale)) + pcauchy(0, 
        location, scale), location, scale)
}
