$(document).ready(function(){
    //Send login information to server
    //Not needed use HTML form action instead
    /*
    $("#signinbutton").click(function(){
        $.post("/login",
            $("#loginForm").serialize(),
            function(response){
                $("#content").html(response)
            });
        return false;
    });
    */

    $("#newuserbutton").click(function(){
        $.post("/createUser",
            $("#loginForm").serialize(),
            function(response){
                $("#content").html(response)
            });
        return false;
    })
});