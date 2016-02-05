$(document).ready(function(){
    //Send login information to server
    $("#loginForm").submit(function(){
        $.post("/login",
            $("#loginForm").serialize(),
            function(response){
                $("#content").html(response)
            });
        return false;
    });
});