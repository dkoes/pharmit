<?php
session_start();

//configuration variable
$db_user = "pharmit";
$db_host = "localhost";
$db_name = "pharmitusers";
$debug = 1;


//this is called for errors that should not be exposed to the user
function fail($pub, $pvt = '')
{
	$msg = $pub;
	if ($debug && $pvt !== '')
		$msg .= ": $pvt";

	//dkoes - try to only fail before generating html header..
	header('Content-Type: text/plain');

	exit("An error occurred ($msg).\n");
}

//produce html with error message
function failhtml($msg)
{
	echo("<html><title>Error</title><body>");
	printf("<h1>%s<h1>", $msg);
	echo("</body></html>");
	exit();
}

//from https://gist.github.com/tylerhall/521810
function generateStrongPassword($length = 9, $add_dashes = false, $available_sets = 'luds')
{
	$sets = array();
	if(strpos($available_sets, 'l') !== false)
		$sets[] = 'abcdefghjkmnpqrstuvwxyz';
	if(strpos($available_sets, 'u') !== false)
		$sets[] = 'ABCDEFGHJKMNPQRSTUVWXYZ';
	if(strpos($available_sets, 'd') !== false)
		$sets[] = '23456789';
	if(strpos($available_sets, 's') !== false)
		$sets[] = '!@#$%&*?';

	$all = '';
	$password = '';
	foreach($sets as $set)
	{
		$password .= $set[array_rand(str_split($set))];
		$all .= $set;
	}

	$all = str_split($all);
	for($i = 0; $i < $length - count($sets); $i++)
		$password .= $all[array_rand($all)];

	$password = str_shuffle($password);

	if(!$add_dashes)
		return $password;

	$dash_len = floor(sqrt($length));
	$dash_str = '';
	while(strlen($password) > $dash_len)
	{
		$dash_str .= substr($password, 0, $dash_len) . '-';
		$password = substr($password, $dash_len);
	}
	$dash_str .= $password;
	return $dash_str;
}

//html boilerplate for regular page
function headerhtml()
{
?>
<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8">
<meta http-equiv="content-script-type" content="text/javascript">
<meta http-equiv="content-style-type" content="text/css">
<link rel="stylesheet" type="text/css" href="create.css" />
<link href='http://fonts.googleapis.com/css?family=Open+Sans' rel='stylesheet' type='text/css'>

<title>Pharmit Library Creation</title>
</head>

<body>

<?php
}

function footerhtml()
{
	echo('</body></html>');
}


//if we aren't logged in already, show login screen
if (!isset($_REQUEST["op"]) && !isset($_SESSION['userid']))
{
	headerhtml();
?>
<div class="loginpage">
<span class="font-3">log in</span><br>
<span class="font">log in to build or manage public or private libraries.</span><br><br>
<div class="loginbox">
<form action="create.php" method="POST">
<input type="hidden" name="op" value="login">
<span class="font-2">email:</span>
<input type="text" autofocus="autofocus" name="user" size="60" class="input"><br><br>
<span class="font-2">password:</span>
<input type="password" name="pass" size="60" class="input"><br><br>
<input type="submit" value="log in" class="submit">
</form>
</div>
<br><br><br>
<span class="font-2">If you don't have an account or have lost your password you can <a href="create.php?op=register">register</a> for one
or <a href="create.php?op=guestlogin">log in as a guest</a>.</span>
<br><br><span class="font-2"><a href="index.php">return to pharmit</a></span><br>

</div>
<?php
	footerhtml();
}
else if(isset($_REQUEST["op"])) //operation
{
	$op = $_REQUEST["op"];
	switch ($op) {
		case "login":
			//check user/pass
			$user = $_POST['user'];
			$pass = $_POST['pass'];

			$db = new mysqli($db_host, $db_user, "", $db_name);
			if (mysqli_connect_errno())
				fail('MySQL connect', mysqli_connect_error());

			($stmt = $db->prepare('SELECT password FROM users WHERE email=?')) ||
				fail('Prepare users', $db->error);
			$stmt->bind_param('s', $user) || fail('Bind user', $db->error);
			$stmt->execute();
			$stmt->store_result();
			if($stmt->num_rows > 0) { //have valid username
				$correctpass = "";
				$stmt->bind_result($correctpass) || fail('Bind pass', $db->error);
				if(!$stmt->fetch() && $db->errno)
					fail('Fetch pass', $db->error);

				if($correctpass == $pass) {
					session_regenerate_id();
					$_SESSION['userid']  = $user;
					session_write_close();
					header("location:create.php");
					exit();

				} else {
					failhtml("Invalid password for $user");
				}

			} else {
				failhtml("Invalid user");
			}
			break;
		case "guestlogin":
			$_SESSION['userid']  = "guest";
			header("location:create.php"); //reload now that session is set with no op
			exit();
			break;
		case "register":

			headerhtml();
			?>
			<div class="loginpage">
			<span class="font-3">register</span><br>
			<span class="font">provide your information and we will email you a password.</span><br><br>
			<div class="loginbox">
			<form action="create.php" method="POST">
			<input type="hidden" name="op" value="doregister">
			<span class="font-2">email:</span>
			<input type="text" autofocus="autofocus" name="email" size="60" class="input"><br><br>
			<span class="font-2">name:</span>
			<input type="text" name="name" size="60" class="input"><br><br>
			<span class="font-2">institution:</span>
			<input type="text" name="place" size="60" class="input"><br><br>
			<input type="submit" value="submit registration" class="submit">
			</form>
			</div>
			<br><br><span class="font-2"><a href="index.php">return to pharmit</a></span><br>

			</div>
			<?php
				footerhtml();
			break;
		case "doregister":
			//create a user
			$email = $_POST['email'];
			$name = $_POST['name'];
			$place = $_POST['place'];

			if (!filter_var($email, FILTER_VALIDATE_EMAIL)) {
				failhtml("Invalid email: $email");
			}
			//generate a password - we assign a password and don't let the
			//user change since this way we (well, I) don't feel bad about
			//storing everything in clear text
			$pass = generateStrongPassword();

			//insert user (or replace, this is what you do when you lose your password)
			$db = new mysqli($db_host, $db_user, "", $db_name);
			if (mysqli_connect_errno())
				fail('MySQL connect', mysqli_connect_error());

			$stmt = $db->prepare('REPLACE INTO users (email, password, name, institution) VALUES (?,?,?,?)');
			$stmt->bind_param('ssss', $email,$pass,$name,$place) || fail('Bind register', $db->error);
			if (!$stmt->execute()) {
					failhtml('Unexpected error registering. Please try again later or contact the site administrator. ');
			}
			else {
				mail( $email, "Pharmit Password" ,
				"Your password is:\n$pass\n\nIf you lose this password you can simply re-register with the same email.");
				echo("Your password has been mailed to $email.  Please check your spam filters.");
			}
			$stmt->close();

			break;
		case "create":
			//have to validate all the inputs
			headerhtml();
			echo("Under construction...");
			footerhtml();
			break;
		case "status":
			//have to validate all the inputs
			headerhtml();
			echo("Under construction...");
			footerhtml();
			break;
		case "logout":
			//remove session totally

			// Unset all of the session variables.
			session_unset();
			// If it's desired to kill the session, also delete the session cookie.
			// Note: This will destroy the session, and not just the session data!
			if (ini_get("session.use_cookies")) {
				$params = session_get_cookie_params();
				setcookie(session_name(), '', time() - 42000,
				$params["path"], $params["domain"],
				$params["secure"], $params["httponly"]
				);
			}

			// Finally, destroy the session.
			session_destroy();
			//back to login screen
			header("location:create.php");

			break;
	}
}
else //logged in, let's create some databases
{
	headerhtml();
	?>
	<div class="createpage">
	<span class="font-3">create</span><br>
	<div class="loginbox">
	<form action="create.php" method="POST">

	<input type="hidden" name="op" value="create">
	<span class="font-2">a short descriptive name of the database:</span><br>
	<input type="text" autofocus="autofocus" name="dbname" size="60" class="input-2"><br><br>
	<span class="font-2">a longer description: </span><br>
	<textarea rows="4" cols="50" name="description">
</textarea><br>
	<span class="font-4">please include any information you think may be useful (supports HTML markup)</i></span><br><br>

	<span class="font-2">access:</span><br>
 	<select name="access">
  		<option value="public">public - anyone will be able to view and search</option>
  		<option value="private">private - a passcode will be required to view and search</option>
	</select><br>
  	<span class="font-4">there are additional limitations on the number and size of private databases.</span><br><br>

	<span class="font-2">compound file (either .smi or .sdf.gz):</span><br>
	<input type="file" name="compounds">
	<br>
	<input type="submit" value="submit" class="submit">
	</form>
	</div>


	</div>
	<br><br>
	<form action="create.php" method="POST">
	<input type="hidden" name="op" value="logout">
	<input type="submit" value="log out" class="submit-2">
	</form>
<?php

	footerhtml();
}
?>

