function Menu(){}

Menu.prototype.loadLoginPage = function(){
    $("#content").load("login.html");
};

Menu.prototype.toggleLoginButton = function(){
    $("#loginButton").toggle();
    $("#logoutButton").toggle();
};

$(document).ready(function(){
    window.menu = new Menu();

    //Listeners
    $("#loginButton").click(function(){
        window.menu.loadLoginPage();
        window.menu.toggleLoginButton();
    });

    $("#logoutButton").click(function(){
        window.location.href="logout";
        window.menu.toggleLoginButton();
    });
});


