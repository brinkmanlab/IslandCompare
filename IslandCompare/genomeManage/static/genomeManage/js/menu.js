function Menu(){}

Menu.prototype.loadLoginPage = function(){
    $("#content").load("login.html");
};

$(document).ready(function(){
    window.menu = new Menu();

    //Listeners
    $("#loginButton").click(function(){
        window.menu.loadLoginPage();
    });

    $("#logoutButton").click(function(){
        window.location.href="logout";
    });
});


