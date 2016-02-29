$(document).ready(function(){
    //loadGenomesToTable(); //replaced with DataTables (jquery plugin)
    //loadJobsToTable(); //replaced with DataTables

    //Setup Listeners below

    //Send Job Request to Server When submit on genome list form is clicked
    $("#genomeListForm").submit(function(){
        //Disabled button to prevent multiple submission of same job
        DisableAnalysisButton();

        //Get Selected Rows in the Table
        var selectedData = $("#genomeTable").DataTable().rows( { selected: true } ).data();
        var runList = [];
        for (var rowIndex=0;rowIndex<selectedData.length;rowIndex++){
            runList.push(selectedData[rowIndex][0]);
        }
        //Serialize the array and add the array of selected genome ids to the array
        var values = $("#genomeListForm").serializeArray();
        values.push({
            name: "selectedSequences",
            value: runList
        });
        values = jQuery.param(values);
        //Send the serialized array to the server
        $.post("/submitJob",
            values,
            function(response){
                ReloadJobsTable();
            });
        return false;
    });

    //Enable Submission of jobs after genome is clicked
    $("#genomeListForm").on('click',"tr",null,function(){
        EnableAnalysisButton();
    })
});

function EnableAnalysisButton(){
    $("#runAnalysisButton").prop("disabled",false);
    $("#runAnalysisButton").removeClass("disabled")
}

function DisableAnalysisButton(){
    $("#runAnalysisButton").prop("disabled",true);
    $("#runAnalysisButton").addClass("disabled")
}

function ReloadGenomesTable(){
    var genomeTable = $("#genomeTable").dataTable();
    $.ajax({
        type:"GET",
        url:"getGenomes",
        success: function(data){
            genomeTable.fnClearTable();
            genomeTable.fnAddData(data['data']);
        }
    })
}

function ReloadJobsTable(){
    var jobsTable = $("#jobTable").dataTable();
    $.ajax({
        type:"GET",
        url:"getJobs",
        success: function(data){
            jobsTable.fnClearTable();
            jobsTable.fnAddData(data['data']);
        }
    })
}

// Load Genomes into Status Table in UI (No Longer Used)
function loadGenomesToTable(){
    $.ajax({
        type:"GET",
        url: "/getGenomes",
        success: function(data){
            genomeTable = $("#genomeTable > tbody");
            $("#genomeTable tbody > tr").remove(); // Clear all data from UI, this will be reloaded
            for (var index=0;index<data.length;index++){
                delete data[index]['uploader']; // Remove uploader from output dict
                tablerowbuilder = "";
                tablerowbuilder = tablerowbuilder.concat("<tr>");
                genome = data[index];
                tablerowbuilder = tablerowbuilder.concat("<td>"+data[index]["name"]+"</td>");
                tablerowbuilder = tablerowbuilder.concat("<td>"+data[index]["length"]+"</td>");
                tablerowbuilder = tablerowbuilder.concat("<td>"+data[index]["description"]+"</td>");
                tablerowbuilder = tablerowbuilder.concat("<td>"+data[index]["uploadedName"]+"</td>");
                tablerowbuilder = tablerowbuilder.concat("<td>"+
                    "<input type=\"checkbox\" name=\"jobCheckList\" value="+data[index].id+" />"+"</td>");
                tablerowbuilder = tablerowbuilder.concat("</tr>");
                genomeTable.append(tablerowbuilder);
            }
        }
    })
}

// Load Jobs into the Job Table in UI (No Longer Used)
function loadJobsToTable(){
    $.ajax({
        type:"GET",
        url:"getJobs",
        success: function(data){
            jobsTable = $("#jobTable > tbody");
            $("#jobTable tbody > tr").remove(); // Clear all data from UI, this will be reloaded
            for (var index=0;index<data.length;index++){
                tablerowbuilder = "";
                tablerowbuilder = tablerowbuilder.concat("<tr>");
                for (var item in data[index]){
                    tablerowbuilder = tablerowbuilder.concat("<td>"+data[index][item]+"</td>");
                }
                tablerowbuilder = tablerowbuilder.concat("<td>"+"<a href=\'"+"getAlignment?id="+
                    data[index].id+"\'"+">Mauve</a></td>");
                tablerowbuilder = tablerowbuilder.concat("</tr>");
                jobsTable.append(tablerowbuilder);
            }
        }
    })
}