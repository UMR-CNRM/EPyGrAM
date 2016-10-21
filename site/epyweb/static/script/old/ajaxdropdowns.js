$(document).ready(function() {
	/*
	 * The jQuery $(document).ready function can be placed on multiple script files.
	 * 
	 * It's smart enough to figure it out, and they'll all be executed 
	 * WHEN the page has loaded.
	 * 
	 */
	
	// this will bind the change event of our countries drop down
	$("#countries").change(function() {
		//alert($(this).val());

		// jquery ajax args must be a JSON string.  So, you can build a string like this:
		// (see the next ajax call for a different approach)
		var args = '{"file":"' + $(this).val() + '"}';

		//this call will get back a JSON list, and use client side script to populate it
		$.ajax({
			type : "POST",
			async : false,
			url : "/getfieldsasjson",
			data : args,
			contentType : "application/json; charset=utf-8",
			dataType : "json",
			success : function(regions) {
				$("#regions").empty();
				$.each(regions, function(index, region) {
					$("#regions").append("<option value='" + region + "'>" + region + "</option>");
				});
			},
			error : function(response) {
				alert(response.responseText);
			}
		});


		// if you'd rather take an object approact for args
		args = {}
		args.country = $(this).val();
		
		//IF you build your args as an object, it must be serialized into json before submission...
		jsargs = JSON.stringify(args);

		//this call will get back HTML that was drawn on the server, and simply inject it into the proper select element
		$.ajax({
			type : "POST",
			async : false,
			url : "/getregionsashtml",
			data : jsargs,
			contentType : "application/json; charset=utf-8",
			dataType : "html",
			success : function(response) {
				$("#regions2").html(response);
			},
			error : function(response) {
				alert(response.responseText);
			}
		});

		//this call gets 100% of it's content all generated server side
		// if you like building stuff in Python instead of javascript.
		$.ajax({
			type : "POST",
			async : false,
			url : "/fastat",
			data : args,
			contentType : "application/json; charset=utf-8",
			dataType : "html",
			success : function(response) {
				$("#coolstuff").html(response);
			},
			error : function(response) {
				alert(response.responseText);
			}
		}); 

	});
});
