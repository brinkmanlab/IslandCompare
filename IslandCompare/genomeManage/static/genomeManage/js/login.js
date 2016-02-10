$(document).ready(function(){
    //Send login information to server
    $("#signinbutton").click(function(){
        $.post("/login",
            $("#loginForm").serialize(),
            function(response){
                $("#content").html(response)
            });
        return false;
    });

    $("#newuserbutton").click(function(){
        $.post("/createUser",
            $("#loginForm").serialize(),
            function(response){
                $("#content").html(response)
            });
        return false;
    })
});